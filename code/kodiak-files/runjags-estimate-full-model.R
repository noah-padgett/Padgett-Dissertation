# Estimate Full Model
# Load packages and sim_data function
source("code/load_packages.R")
# packages needed for this specific file
# mvtnorm
# runjags

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
  C = paravec[3],
  tune = 50
)

jags.model <- "model{

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
      #mu.rt[p,i]<- icept[i] - speed[p] - rho * abs(dev[p,i])
      lrt[p,i] ~ dnorm(icept[i] - speed[p] - rho * abs(dev[p,i]), prec[i])

      # compute ELRT
      logit(ELRT[p,i]) = exp(icept[i] - speed[p])/abs(dev[p,i])
      # versusi
      #logit(ELRT[p,i]) = exp(icept[i] - speed[p])*abs(dev[p,i])

      # MISCLASSIFICATION MODEL
      for(c in 1:C){
        # generate informative prior for misclassificaiton
        #   parameters based on RT
        for(ct in 1:C){
          alpha[p,i,ct,c] <- ifelse(c == ct,
                                    ELRT[p,i]*tune,
                                    (1/(C-1))*(1-ELRT[p,i])*tune
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

jags.params.s <- c(
  # ability measurement model
  "tau", "lambda",
  # speed measurement parameters
  "rho",  "icept","prec.s", "sigma.ts"
)

fit.runjags <- run.jags(
  model=jags.model
  , data=mydata
  , keep.jags.files=F
  #, tempdir = '/data/padgettn/temp/'
  , monitor = jags.params.s
  , n.chains = 4
  , method = "rjags"
  , sample = 5000
  , adapt= 2000
)

summary(fit.runjags)



