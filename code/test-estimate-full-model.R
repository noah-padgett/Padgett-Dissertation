# Estimate Full Model
# Load packages and sim_data function
source("code/load_packages.R")
source("code/simulate_data.R")


# get data
paravec <- c(
  N = 500, J = 5, C = 4,
  etaCor = .23, etasd1 = 1, etasd2 = sqrt(0.1),
  lambda=0.9, nu=1.5, sigma.ei=0.25, rho1=0.1)
sTau <- matrix(
  c(-2.5, -0.5, 1,
    -2, -0.5, 1.5,
    -2.25, -0.5, 1.5,
    -1.75, -0.5, 2,
    -1.75, -0.5, 1.75),
  ncol=paravec[3]-1, nrow=paravec[2], byrow=T
)
sim.data <- simulate_data_misclass(paravec, tau=sTau)

mydata <- list(
  y = sim.data$y,
  lrt = sim.data$logt,
  N = nrow(sim.data$y),
  nit=ncol(sim.data$y),
  C = paravec[3]
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dcat(omega[p,i, ])

      # LRV
      ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)

      # Pr(nu = 3)
      pi[p,i,5] = phi(ystar[p,i] - tau[i,4])
      # Pr(nu = 4)
      pi[p,i,4] = phi(ystar[p,i] - tau[i,3]) - phi(ystar[p,i] - tau[i,4])
      # Pr(nu = 3)
      pi[p,i,3] = phi(ystar[p,i] - tau[i,2]) - phi(ystar[p,i] - tau[i,3])
      # Pr(nu = 2)
      pi[p,i,2] = phi(ystar[p,i] - tau[i,1]) - phi(ystar[p,i] - tau[i,2])
      # Pr(nu = 1)
      pi[p,i,1] = 1 - phi(ystar[p,i] - tau[i,1])

      # log-RT model
      dev[p,i]<-lambda[i]*(eta[p] - (tau[i,1]+tau[i,2]+tau[i,3]+tau[i,4])/4)
      mu.rt[p,i]<- icept[i] - speed[p] - rho * abs(dev[p,i])
      lrt[p,i] ~ dnorm(mu.rt[p,i], prec[i])

      # compute ELRT
      logit(ELRT[p,i]) = exp(icept[i] - speed[p])/abs(dev[p,i])

      # MISCLASSIFICATION MODEL
      for(c in 1:C){
        # generate informative prior for misclassificaiton
        #   parameters based on RT
        for(ct in 1:C){
          alpha[p,i,ct,c] <- ifelse(c == ct,
                                    ELRT[p,i]*10,
                                    (1/(C-1))*(1-ELRT[p,i])*10
          )
        }
        # sample misclassification parameters using the informative priors
        gamma[p,i,c,1:C] ~ ddirch(alpha[p,i,c,1:C])
        # observed category prob (Pr(y=c))
        omega[p,i, c] = gamma[p,i,c,1]*pi[p,i,1] +
          gamma[p,i,c,2]*pi[p,i,2] +
          gamma[p,i,c,3]*pi[p,i,3] +
          gamma[p,i,c,4]*pi[p,i,4] +
          gamma[p,i,c,5]*pi[p,i,5]
      }

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
    speed[p]~dnorm(sigma.ts*eta[p],prec.s)  # latent speed
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s ~ dgamma(.1,.1)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # Thresholds
    tau[i, 1] ~ dnorm(0.0,0.1)
    tau[i, 2] ~ dnorm(0, 0.1);T(tau[i, 1],)
    tau[i, 3] ~ dnorm(0, 0.1);T(tau[i, 2],)
    tau[i, 4] ~ dnorm(0, 0.1);T(tau[i, 3],)
    # loadings
    lambda[i] ~ dnorm(0, .44);T(0,)
    # LRV total variance
    # total variance = residual variance + fact. Var.
    theta[i] = 1 + pow(lambda[i],2)
  }
  rho~dnorm(0,.1);I(0,)

  # important parameters
  sigma.t = pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts = sigma.ts/(pow(sigma.t,0.5))

  # omega
  l2 = (lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]+lambda[6]+lambda[7]+lambda[8]+lambda[9]+lambda[10]+lambda[11]+lambda[12]+lambda[13]+lambda[14]+lambda[15]+lambda[16]+lambda[17]+lambda[18]+lambda[19]+lambda[20])
  reli.omega = (pow(l2,2))/(pow(l2,2)+ nit)
}
jags.params <- c(
  # ability measurement model
  "tau", "lambda", "theta", "pi",
  # speed measurement parameters
  "rho",  "icept","prec.s", "sigma.t", "sigma.ts", "cor.ts",
  # misclassification parameters
  "omega", "gamma", "ELRT",
  # reliability
  "reli.omega"
)


# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = jags.params,
  n.iter=2
  , n.chains = 1
  #,n.thin = 50
  #,n.burnin = 0000
  #,n.chains = 2
)



# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
fit.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit.mcmc) <- c("chain", "iter", a)


mcmc_areas(fit.mcmc,pars = "reli.omega", prob = 0.8)
mcmc_areas(fit.mcmc,regex_pars = "lambda", prob = 0.8)

# get a prior predicted distribution of reli
y.null <- lrt.null <- matrix(NA, ncol=5, nrow=200)
mydata <- list(
  y = y.null,
  lrt = lrt.null,
  N   = 200,
  nit = 5,
  C   = 3
)
# fit model
fit.0 <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = jags.params,
  n.iter=2000,
  #n.thin = 50,
  #n.burnin = 0000,
  n.chains = 2
)

# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit.0)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
fit.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit.mcmc) <- c("chain", "iter", a)



# prior predictive distribution of reli.
niter <- 1000
ppd.omega <- numeric(niter)
for(iter in 1:niter){
  nit <- 5
  lambda <- theta <- rep(0, nit)
  lambda <- rnorm(nit, 0, 1.5)

  ppd.omega[iter] <- (sum(lambda)**2)/(sum(lambda)**2 + nit)
}

plot(density(ppd.omega))




# test the full model using runjags
paravec <- c(N = 500, J = 10, C = 5, etaCor = .50, etasd1 = 1, etasd2 = sqrt(0.1), lambda=0.9, nu=1.5, sigma.ei=0.25, rho1=0.3)
sim.data <- simulate_data_misclass(paravec)


mydata <- list(
  y = sim.data$y,
  lrt = sim.data$logt,
  N = nrow(sim.data$y),
  nit=ncol(sim.data$y),
  C = 5
)

jags.model <- "model{

  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dcat(omega[p,i, ])

      # LRV
      #ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)

      # Pr(nu = 3)
      pi[p,i,5] = phi(lambda[i]*eta[p] - tau[i,4])
      # Pr(nu = 4)
      pi[p,i,4] = phi(lambda[i]*eta[p] - tau[i,3]) - phi(lambda[i]*eta[p] - tau[i,4])
      # Pr(nu = 3)
      pi[p,i,3] = phi(lambda[i]*eta[p] - tau[i,2]) - phi(lambda[i]*eta[p] - tau[i,3])
      # Pr(nu = 2)
      pi[p,i,2] = phi(lambda[i]*eta[p] - tau[i,1]) - phi(lambda[i]*eta[p] - tau[i,2])
      # Pr(nu = 1)
      pi[p,i,1] = 1 - phi(lambda[i]*eta[p] - tau[i,1])

      # log-RT model
      dev[p,i]<-lambda[i]*(eta[p] - (tau[i,1]+tau[i,2]+tau[i,3]+tau[i,4])/4)
      #mu.rt[p,i]<- icept[i] - speed[p] - rho * abs(dev[p,i])
      lrt[p,i] ~ dnorm(icept[i] - speed[p] - rho * abs(dev[p,i]), prec[i])

      # compute ELRT
      logit(ELRT[p,i]) = exp(icept[i] - speed[p])/abs(dev[p,i])

      # MISCLASSIFICATION MODEL
      for(c in 1:C){
        # generate informative prior for misclassificaiton
        #   parameters based on RT
        for(ct in 1:C){
          alpha[p,i,ct,c] <- ifelse(c == ct,
                                    ELRT[p,i]*50,
                                    (1/(C-1))*(1-ELRT[p,i])*50
          )
        }
        # sample misclassification parameters using the informative priors
        gamma[p,i,c,1:C] ~ ddirch(alpha[p,i,c,1:C])
        # observed category prob (Pr(y=c))
        omega[p,i, c] = gamma[p,i,c,1]*pi[p,i,1] +
          gamma[p,i,c,2]*pi[p,i,2] +
          gamma[p,i,c,3]*pi[p,i,3] +
          gamma[p,i,c,4]*pi[p,i,4] +
          gamma[p,i,c,5]*pi[p,i,5]
      }

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
    speed[p]~dnorm(sigma.ts*eta[p],prec.s)  # latent speed
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s ~ dgamma(.1,.1)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # Thresholds
    tau[i, 1] ~ dnorm(0.0,0.1)
    tau[i, 2] ~ dnorm(0, 0.1)T(tau[i, 1],)
    tau[i, 3] ~ dnorm(0, 0.1)T(tau[i, 2],)
    tau[i, 4] ~ dnorm(0, 0.1)T(tau[i, 3],)
    # loadings
    lambda[i] ~ dnorm(0, .44)T(0,)
  }
  rho~dnorm(0,.1)I(0,)

}"

jags.params.l <- c(
  "pi", "omega",    "gamma"  ,  "ELRT"
)
jags.params.s <- c(
  # ability measurement model
  "tau", "lambda",
  # speed measurement parameters
  "rho",  "icept","prec.s", "sigma.ts"
)
# posterior mutated functions
my.funcs <- function(x){
  # initialize values
  theta1 <- 1 + x[,"lambda[1]"]**2
  theta2 <- 1 + x[,"lambda[2]"]**2
  theta3 <- 1 + x[,"lambda[3]"]**2
  theta4 <- 1 + x[,"lambda[4]"]**2
  theta5 <- 1 + x[,"lambda[5]"]**2

  sigma.t = pow(x[,"prec.s"], -1) + pow(x[,"sigma.ts"], 2)
  cor.ts = x[,"sigma.ts"]/(pow(sigma.t,0.5))

  # omega
  l2 = (x[,"lambda[1]"] + x[,"lambda[2]"] + x[,"lambda[3]"]+ x[,"lambda[4]"] + x[,"lambda[5]"])
  reli.est = (pow(l2,2))/(pow(l2,2)+ 5)
  return(list(reli.est=reli.est, cor.ts=cor.ts,
              theta1=theta1, theta2=theta2, theta3=theta3,
              theta4=theta4,  theta5=theta5))
}

# initial values
my.inits <- function(chain){
  .RNG.seed <- c(1,2)[chain]
  .RNG.name <- c("base::Super-Duper", "base::Wichmann-Hill")[chain]
  return(list(.RNG.seed = .RNG.seed, .RNG.name=.RNG.name))
}
library(runjags)

Sys.time()

fit.runjags <- run.jags(
  model=jags.model
  , data=mydata
  , inits = my.inits
  #, max.time ="30m"
  , keep.jags.files=T
  , monitor = jags.params.s
  #, noread.monitor = jags.params.l
  , n.chains = 2
  , method = "parallel"
  #, mutate = my.funcs
  , sample = 2000
  , adapt= 2000
)
Sys.time()

# I=5 => 2.5hr
# I=10 => 4.5hr

fit.runjags
summary(fit.runjags)






jags.model <- "model{

  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dcat(pi[p,i, ])

      # LRV
      ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)

      # Pr(nu = 5)
      pi[p,i,5] = phi(ystar[p,i] - tau[i,4])
      # Pr(nu = 4)
      pi[p,i,4] = phi(ystar[p,i] - tau[i,3]) - phi(ystar[p,i] - tau[i,4])
      # Pr(nu = 3)
      pi[p,i,3] = phi(ystar[p,i] - tau[i,2]) - phi(ystar[p,i] - tau[i,3])
      # Pr(nu = 2)
      pi[p,i,2] = phi(ystar[p,i] - tau[i,1]) - phi(ystar[p,i] - tau[i,2])
      # Pr(nu = 1)
      pi[p,i,1] = 1 - phi(ystar[p,i] - tau[i,1])

      # log-RT model
      lrt[p,i] ~ dnorm(icept[i] - speed[p] - rho *(lambda[i]*(eta[p] - (tau[i,1]+tau[i,2]+tau[i,3]+tau[i,4])/4)), prec[i])

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
    speed[p]~dnorm(sigma.ts*eta[p],prec.s)  # latent speed
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s ~ dgamma(.1,.1)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # Thresholds
    tau[i, 1] ~ dnorm(0.0,0.1)
    tau[i, 2] ~ dnorm(0, 0.1)T(tau[i, 1],)
    tau[i, 3] ~ dnorm(0, 0.1)T(tau[i, 2],)
    tau[i, 4] ~ dnorm(0, 0.1)T(tau[i, 3],)
    # loadings
    lambda[i] ~ dnorm(0, .44)T(0,)
  }
  rho~dnorm(0,.1)I(0,)

}"

jags.params <- c(
  # ability measurement model
  "tau", "lambda", "pi",
  # speed measurement parameters
  "rho",  "icept","prec.s", "sigma.t", "sigma.ts"
)

# initial values
my.inits <- function(){
  list(
    "eta"=rnorm(data$N, 0, 1)
  )
}

fit.runjags <- autorun.jags(
  model=jags.model
  , data=mydata
  , inits = my.inits
  , max.time ="2h"
  , monitor=jags.params
  , n.chains = 4
  , method = "parallel"
)

fit.runjags





# estimate extroversion model
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
    Mt = mean(Time, na.rm=T)
  )

mydata <- list(
  y = extraversion[,1:10],
  lrt = log(extraversion[,11:20]),
  N = nrow(extraversion),
  nit=10,
  C=2
)

jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dbern(omega[p,i,2])

      # LRV
      ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)

      # Pr(nu = 1)
      pi[p,i,2] = phi(ystar[p,i] - tau[i])
      pi[p,i,1] = 1 - phi(ystar[p,i] - tau[i])

      # log-RT model
      dev[p,i]<- lambda[i]*(eta[p] - tau[i])
      mu.rt[p,i]<- icept[i] - speed[p] - rho * abs(dev[p,i])
      lrt[p,i] ~ dnorm(mu.rt[p,i], prec[i])

      # compute ELRT
      ELRT[p,i] = exp(icept[i] - speed[p])

      # MISCLASSIFICATION MODEL
      for(c in 1:C){
        # generate informative prior for misclassificaiton
        #   parameters based on RT
        for(ct in 1:C){
          alpha[p,i,ct,c] <- ifelse(c == ct,
                                    ilogit(ELRT[p,i])*10,
                                    (1/(C-1))*(1-ilogit(ELRT[p,i]))*10
          )
        }
        # sample misclassification parameters using the informative priors
        gamma[p,i,c,1:C] ~ ddirch(alpha[p,i,c,1:C])
        # observed category prob (Pr(y=c))
        omega[p,i, c] = gamma[p,i,c,1]*pi[p,i,1] +
          gamma[p,i,c,2]*pi[p,i,2]
      }

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
    speed[p]~dnorm(sigma.ts*eta[p],prec.s)  # latent speed
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s ~ dgamma(.1,.1)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # Thresholds
    tau[i] ~ dnorm(0.0,0.1)
    # loadings
    lambda[i] ~ dnorm(0, .44);T(0,)
    # LRV total variance
    # total variance = residual variance + fact. Var.
    theta[i] = 1 + pow(lambda[i],2)
  }
  negrho =- -1*rho
  rho~dnorm(0,.1);I(0,)

  # important parameters
  sigma.t = pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts = sigma.ts/(pow(sigma.t,0.5))

  # omega
  l2 = (lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]+lambda[6]+lambda[7]+lambda[8]+lambda[9]+lambda[10])
  reli.est = (pow(l2,2))/(pow(l2,2)+ nit)
}
jags.params <- c(
  # ability measurement model
  "tau", "lambda", "theta", "pi",
  # speed measurement parameters
  "rho",  "icept","prec.s", "sigma.t", "sigma.ts", "cor.ts",
  # misclassification parameters
  "omega", "gamma", "ELRT",
  # reliability
  "reli.est"
)


# fit model
fit <- jags(
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=2000,
  n.chain=2,
  model.file=jags.model)

fit$BUGSoutput$summary[
  !(rownames(fit$BUGSoutput$summary) %like% "ELRT" |
      rownames(fit$BUGSoutput$summary) %like% "gamma"|
      rownames(fit$BUGSoutput$summary) %like% "pi" |
      rownames(fit$BUGSoutput$summary) %like% "omega"),]

library(ltm)
?ltm::tpm()
fit.3pl <- ltm::tpm(extraversion[,1:10])
fit.3pl

fit.2pl <- ltm::tpm(extraversion[,1:10], max.guessing = 0)
fit.2pl






# trying to match LTM
mydata <- list(
  y = extraversion[,1:10],
  N = nrow(extraversion),
  nit=10
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dbern(pi[p,i])

      # LRV
      ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)

      pi[p,i] = phi(ystar[p,i] - tau[i])

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
  }
  for(i in 1:nit){
    # Thresholds
    tau[i] ~ dnorm(0.0,0.1)
    tau.d[i] = tau[i]/lambda[i] # rescale tau
    # loadings
    lambda[i] ~ dnorm(0, .44);T(0,)
    # LRV total variance
    # total variance = residual variance + fact. Var.
    theta[i] = 1 + pow(lambda[i],2)
  }
}
jags.params <- c(
  # ability measurement model
  "tau","tau.d", "lambda", "theta"
)


# fit model
fit <- jags(
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=2000,
  n.chain=2,
  model.file=jags.model)
print(fit)
fit.auto <- autojags(fit, n.update = 10, n.iter = )
print(fit)


library(runjags)

jags.model <- "model{
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dbern(pi[p,i])

      # LRV
      ystar[p,i] ~ dnorm(lambda[i]*eta[p], 1)

      pi[p,i] = phi(ystar[p,i] - tau[i])

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
  }
  for(i in 1:nit){
    # Thresholds
    tau[i] ~ dnorm(0.0,0.1)
    tau.d[i] = tau[i]/lambda[i] # rescale tau
    # loadings
    lambda[i] ~ dnorm(0, .44)T(0,)
    # LRV total variance
    # total variance = residual variance + fact. Var.
    theta[i] = 1 + pow(lambda[i],2)
  }
}
"
fit.runjags <- autorun.jags(
  model=jags.model
  , data=mydata
  , max.time ="5m"
  , monitor=jags.params
  , n.chains = 4
  , method = "parallel"
)

fit.runjags
