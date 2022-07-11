# test with extraversion data - diffIRT package
source("code/load_packages.R")
library(diffIRT)
options(max.print = 10000, scipen = 10, digits=2)

# useful functions
invlogit <- function(x) {exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

data("extraversion")

d1 <- extraversion %>%
  as.data.frame() %>%
  select(contains("X"))%>%
  mutate(id = 1:n()) %>%
  pivot_longer(
    cols=contains("X"),
    names_to = c("item"),
    values_to = "Response"
  ) %>%
  mutate(
    item = ifelse(nchar(item)==4,substr(item, 3,3),substr(item, 3,4))
  )
d2 <- extraversion %>%
  as.data.frame() %>%
  select(contains("T"))%>%
  mutate(id = 1:n()) %>%
  pivot_longer(
    cols=starts_with("T"),
    names_to = c("item"),
    values_to = "Time"
  ) %>%
  mutate(
    item = ifelse(nchar(item)==4,substr(item, 3,3),substr(item, 3,4))
  )
dat <- left_join(d1, d2)

dat %>%
  select(item, Response, Time) %>%
  group_by(item) %>%
  summarize(
    M1 = mean(Response, na.rm=T),
    Mt = mean(Time, na.rm=T),
    SDt = sd(Time, na.rm=T)
  )
# eRm
library(eRm)
fit.R <- eRm::RM(extraversion[,1:10])
summary(fit.R)
pp <- person.parameter(fit.R)
itemfit(pp)

#sirt
library(sirt)
mod2 <- sirt::mcmc.2pno(extraversion[,1:10])
#mod2
summary(mod2)

# mirt - 4pl model
library(mirt)
fit <- mirt(extraversion[,1:10], 1, itemtype = rep("4PL",10))
coef(fit)

# Rasch model
mydata <- list(
  x = extraversion[,1:10],
  N = nrow(extraversion),
  nit=10
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # Rasch model
      logit(Prob[p,i]) <- theta[p]-b[i]
      # observable dist.
      x[p,i] ~ dbern(Prob[p,i])
    }
  }
  for(p in 1:N){
    theta[p]~dlogis(0, 1)
  }
  #mu ~ dnorm(0, 0.1)

  b[1] <- -sum(b[2:nit])
  for(i in 2:nit){
    # location parameters
    b[i] ~ dnorm(0, .1)
  }
  # Rasch item fit statistics
  # deviations
  for(i in 1:nit){
    for(p in 1:N){
      w[p,i] <- Prob[p,i]*(1-Prob[p,i])
      z[p,i] <- (x[p,i]-Prob[p,i])/pow(w[p,i], 0.5)
      wz2[p,i] <- w[p,i]*pow(z[p,i],2)
    }
    outfitMSQ[i] <- sum(pow(z[,i], 2))/N
    infitMSQ[i] <- sum(wz2[,i])/sum(w[,i])
  }

} # end jags model


# params to save
jags.params <- c("b", "outfitMSQ", "infitMSQ")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = jags.params,
  n.iter=10000,
  #n.thin = 50,
  #n.burnin = 0000,
  n.chains = 4
)
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(plot.data, regex_pars = "b", prob = 0.8)
mcmc_areas(plot.data, regex_pars = "outfit", prob = 0.8)
mcmc_areas(plot.data, regex_pars = "infit", prob = 0.8)

# set up Rasch response time model
mydata <- list(
  x = extraversion[,1:10],
  lrt = log(extraversion[,11:20]),
  N = nrow(extraversion),
  nit=10
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # Rasch model
      logit(Prob[p,i]) <- theta[p]-b[i]
      # observable dist.
      x[p,i] ~ dbern(Prob[p,i])

      # log-RT model
      dev[p,i]<-theta[p] - b[i]
      mu.rt[p,i]<- icept[i] - speed[p] - abs(dev[p,i])
      lrt[p,i]~dnorm(mu.rt[p,i],prec[i])
    }
  }
  for(p in 1:N){
    theta[p]~dlogis(0,1)
    speed[p]~dnorm(sigma.ts*theta[p],prec.s)
  }

  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # location parameters
    b[i] ~ dnorm(0, .5)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)
  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))
  # Rasch item fit statistics
  # deviations
  for(i in 1:nit){
    for(p in 1:N){
      w[p,i] <- Prob[p,i]*(1-Prob[p,i])
      z[p,i] <- (x[p,i]-Prob[p,i])/pow(w[p,i], 0.5)
      wz2[p,i] <- w[p,i]*pow(z[p,i],2)
    }
    outfitMSQ[i] <- sum(pow(z[,i], 2))/N
    infitMSQ[i] <- sum(wz2[,i])/sum(w[,i])
  }

} # end jags model


# params to save
jags.params <- c("b", "icept","prec","prec.s", "sigma.t", "sigma.ts", "cor.ts", "outfitMSQ", "infitMSQ")


# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = jags.params,
  n.iter=10000,
  #n.thin = 50,
  #n.burnin = 0000,
  n.chains = 4
)
fit2 <- fit
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(
  plot.data,
  pars = c(paste0("b[",1:10,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("icept[",1:10,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("prec[",1:10,"]")),
  prob = 0.8) +
  plot_title


mcmc_areas(
  plot.data,
  pars = c("sigma.t", "sigma.ts", "cor.ts"),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("outfitMSQ[",1:10,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("infitMSQ[",1:10,"]")),
  prob = 0.8) +
  plot_title



# misclassification model
extra <- na.omit(extraversion)
mydata <- list(
  x = extra[,1:10],
  lrt = log(extra[,11:20]),
  N = nrow(extra),
  nit=10, ncat=2
)

jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # observable dist.
      x[p,i] ~ dbern(omega[p,i,2])
      # latent response dist.
      xstar[p,i] ~ dnorm(lambda[i]*eta[p], 1)
      # probs
      Prob[p,i,1] <- phi(xstar[p,i]-tau[i])
      Prob[p,i,2] <- 1 - Prob[p,i,1]
      # log-RT model
      dev[p,i]<- lambda[i]*(eta[p] - tau[i])
      mu.rt[p,i]<- icept[i] - speed[p] - rho*abs(dev[p,i])
      lrt[p,i]~dnorm(mu.rt[p,i],prec[i])
      # ELRT
      #ELRT[p,i] <- ifelse(abs(dev[p,i]) < 0.01, 5, (icept[i] - speed[p])/abs(dev[p,i]))
      # misclassificaiton rates
      for(c in 1:ncat){
        for(ct in 1:ncat){
          gamma[p,i,c,ct] <- ifelse(c == ct,
                                    ilogit(lrt[p,i]),
                                    (1/(ncat-1))*(1-ilogit(lrt[p,i]))
          )
        }
      }

      # misclassification weighting
      for(c in 1:ncat){
        omega[p,i, c] = gamma[p,i,1,c]*Prob[p,i,1] +
          gamma[p,i,2,c]*Prob[p,i,2]
      }

    }
  }
  # person parameters
  for(p in 1:N){
    eta[p]~dnorm(0,1)
    speed[p]~dnorm(sigma.ts*eta[p],prec.s)
  }

  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # location parameters
    tau[i] ~ dnorm(0, .5)
    # factor loadings
    lambda[i] ~ dnorm(0, 0.1);T(0,)
    lambda.std[i] = lambda[i]/pow(theta[i],0.5)
    # total latent response variance
    theta[i] = 1 + pow(lambda[i],2)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)
  rho ~ dnorm(0, 0.5);T(0,)
  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))

  # compute omega
  num = lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]+
    lambda[6]+lambda[7]+lambda[8]+lambda[9]+lambda[10]
  omega.r = pow(num,2)/(pow(num, 2)+ (10))
} # end jags model



# vector of all parameters to save
jags.params <- c("lambda", "lambda.std", "tau","eta", "theta", "icept","prec","prec.s", "sigma.t", "sigma.ts", "cor.ts", "omega.r", "rho")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=2000,
  n.chains = 2
)
fit3 <- fit
print(fit)

fit$BUGSoutput$summary[
  !(rownames(fit$BUGSoutput$summary) %like% "eta"),
]

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")

mcmc_areas(plot.data, pars="omega.r")
mcmc_areas(plot.data, pars="rho")
mcmc_areas(plot.data, pars="cor.ts")


mcmc_areas(
  plot.data,
  pars = c(paste0("b[",1:10,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("icept[",1:10,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("prec[",1:10,"]")),
  prob = 0.8) +
  plot_title


mcmc_areas(
  plot.data,
  pars = c("sigma.t", "sigma.ts", "cor.ts"),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("outfitMSQ[",1:10,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("infitMSQ[",1:10,"]")),
  prob = 0.8) +
  plot_title



# comparing results
print(fit1)
print(fit2)
print(fit3)


fit1$BUGSoutput$summary[paste0("outfitMSQ[",1:10,"]"),]
fit2$BUGSoutput$summary[paste0("outfitMSQ[",1:10,"]"),]
fit3$BUGSoutput$summary[paste0("outfitMSQ[",1:10,"]"),]

fit1$BUGSoutput$summary[paste0("infitMSQ[",1:10,"]"),]
fit2$BUGSoutput$summary[paste0("infitMSQ[",1:10,"]"),]
fit3$BUGSoutput$summary[paste0("infitMSQ[",1:10,"]"),]

# plotting for item 1
plot.dat <- data.frame(
  m1 = fit1$BUGSoutput$sims.matrix[,"outfitMSQ[1]"],
  m2 = fit2$BUGSoutput$sims.matrix[,"outfitMSQ[1]"],
  m3 = fit3$BUGSoutput$sims.matrix[,"outfitMSQ[1]"]
) %>%
  pivot_longer(
    cols=everything(),
    names_to = "model",
    values_to = "outfitMSQ"
  )
ggplot(plot.dat, aes(x=outfitMSQ, group=model, color=model))+
  geom_density()+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  geom_vline(xintercept = 1.29, color="red", linetype="dotted")

plot.dat <- data.frame(
  m1 = fit1$BUGSoutput$sims.matrix[,"infitMSQ[1]"],
  m2 = fit2$BUGSoutput$sims.matrix[,"infitMSQ[1]"],
  m3 = fit3$BUGSoutput$sims.matrix[,"infitMSQ[1]"]
) %>%
  pivot_longer(
    cols=everything(),
    names_to = "model",
    values_to = "infitMSQ"
  )
ggplot(plot.dat, aes(x=infitMSQ, group=model, color=model))+
  geom_density()+
  geom_vline(xintercept = 1, color="black", linetype="dashed")+
  geom_vline(xintercept = 1.2, color="red", linetype="dotted")











# Gelman & Hill version of LRV Parameterization

# Rasch model
mydata <- list(
  x = extraversion[,1:10],
  N = nrow(extraversion),
  nit=10
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      z.lo[p,i] <- -100*equals(x[p,i],0)
      z.hi[p,i] <- 100*equals(x[p,i],1)
      z[p,i] ~ dlogis(Xbeta[p,i], 1); T(z.lo[p,i], z.hi[p,i])
      x[p,i] ~ dbin(p.bound[p,i], 1)
      p.bound[p,i] <- max(0, min(1, prob[p,i]))
      logit(prob[p,i]) <- Xbeta[p,i]
      Xbeta[p,i] <- theta[p]-b[i]
    }
  }
  for(p in 1:N){
    theta[p]~dlogis(0, 1)
  }
  #mu ~ dnorm(0, 0.1)

  b[1] <- -sum(b[2:nit])
  for(i in 2:nit){
    # location parameters
    b[i] ~ dnorm(0, .1)
  }
} # end jags model

# params to save
jags.params <- c("b",  "theta")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=10000,
  n.chains = 4
)
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data, regex_pars = "b", prob = 0.8)
mcmc_areas(plot.data, pars = c("theta[1]", "theta[100]"), prob = 0.8)

