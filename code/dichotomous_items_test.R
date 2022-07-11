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


# misclassification model
extra <- na.omit(extraversion)
mydata <- list(
  y = extra[,1:10]
  , lrt = log(extra[,11:20])
  , N = nrow(extra)
  , nit=10
  , ncat=2
  , tune=50
)

jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # observable dist.
      y[p,i] ~ dbern(omega[p,i,2])
      # latent response dist.
      ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)
      # probs
      Prob[p,i,1] = phi(ystar[p,i]-tau[i])
      Prob[p,i,2] = 1 - Prob[p,i,1]
      # log-RT model
      dev[p,i] = lambda[i]*(eta[p] - tau[i])
      lrt[p,i] ~ dnorm(icept[i] - speed[p] - rho*abs(dev[p,i]), prec[i])
      # ELRT
      logit(ELRT[p,i]) = exp(icept[i] - speed[p])*abs(dev[p,i])

      # misclassificaiton priors
      for(c in 1:ncat){
        for(ct in 1:ncat){
          alpha[p,i,c,ct] <- ifelse(
            c == ct,
            ELRT[p,i]*tune,
            (1-ELRT[p,i])*tune
          )

        }
      }

      # misclassification weighting
      for(c in 1:ncat){
        # sample misclassification parameters using the informative priors
        gamma[p,i,c,1:ncat] ~ ddirch(alpha[p,i,c,1:ncat])

        omega[p,i, c] = gamma[p,i,c,1]*Prob[p,i,1] +
          gamma[p,i,c,2]*Prob[p,i,2]
      }

    }
  }
  # person parameters
  for(p in 1:N){
    eta[p]~dnorm(0,1)
    speed[p]~dnorm(sigma.ts*eta[p],prec.s)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)

  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # location parameters
    tau[i] ~ dnorm(0.0,0.1)
    # factor loadings
    lambda[i] ~ dnorm(0, .44);T(0,)
    lambda.std[i] = lambda[i]/pow(theta[i],0.5)
    # total latent response variance
    theta[i] = 1 + pow(lambda[i],2)
  }

  rho ~ dnorm(0, 0.1);T(0,)

  # compute omega
  lambda_sum[1] = lambda[1]
  for(i in 2:nit){
    #lambda_sum (sum factor loadings)
    lambda_sum[i] = lambda_sum[i-1]+lambda[i]
  }
  omega.r = pow(lambda_sum[nit],2)/(pow(lambda_sum[nit], 2)+ (nit))
} # end jags model



# vector of all parameters to save
jags.params <- c("lambda", "lambda.std", "tau","eta", "theta", "icept","prec","prec.s", "sigma.t", "sigma.ts", "cor.ts", "omega.r", "rho", "omega", "gamma", "ELRT", "Prob", "alpha","lambda_sum")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=500,
  n.chains = 2
)
fit3 <- fit
print(fit)

fit$BUGSoutput$summary[paste0("omega.r"),]
fit$BUGSoutput$summary[
  !(rownames(fit$BUGSoutput$summary) %like% "eta" |
      rownames(fit$BUGSoutput$summary) %like% "Prob" |
      rownames(fit$BUGSoutput$summary) %like% "ELRT" |
      rownames(fit$BUGSoutput$summary) %like% "alpha" |
      rownames(fit$BUGSoutput$summary) %like% "gamma" |
      rownames(fit$BUGSoutput$summary) %like% "omega"),
]


fit$BUGSoutput$summary[paste0("omega[1,",1:10,",1]"),]
fit$BUGSoutput$summary[paste0("omega[1,",1:10,",2]"),]
dat[dat$id==1,]


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


