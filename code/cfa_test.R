
library(readxl)
mydata <- read_excel("data/pools/POOLS_data_2020-11-16.xlsx")

use.var <- c(paste0("Q4_",c(3:5,9,11,15,18)))

# model code
jags.model <- function(){

  ########################################
  # Specify the factor analysis measurement model for the observables
  ##############################################
  for (i in 1:n){
    for(j in 1:J){
      x[i,j] ~ dnorm(tau[j] + ksi[i]*lambda[j], (pow(ksi[i],2)+0.0001)*inv.psi[j])    # distribution for each observable
    }
  }


  ##################################
  # Specify the (prior) distribution for the latent variables
  ####################################
  for (i in 1:n){
    ksi[i] ~ dnorm(kappa, inv.phi)  # distribution for the latent variables
  }


  ######################################
  # Specify the prior distribution for the parameters that govern the latent variables
  ###################################
  kappa <- 0              # Mean of factor 1
  inv.phi ~ dgamma(5, 10) # Precision of factor 1
  phi <- 1/inv.phi        # Variance of factor 1


  ########################################
  # Specify the prior distribution for the measurement model parameters
  ########################################
  for(j in 1:J){
    tau[j] ~ dnorm(3, .1)        # Intercepts for observables
    inv.psi[j] ~ dgamma(5, 10) # Precisions for observables
    psi[j] <- 1/inv.psi[j]   # Variances for observables
  }

  lambda[1] <- 1.0              # loading fixed to 1.0
  lambda.std[1] <- lambda[1]*(pow(psi[1]/phi,0.5) )
  for (j in 2:J){
    lambda[j] ~ dnorm(1, .1)    # prior distribution for the remaining loadings
    lambda.std[j] = lambda[j]*(pow(psi[j]/phi,0.5) )# standardized lambda
  }

}

dat <- na.omit(mydata[,use.var])
jags.dat <- list(
  n = nrow(dat), J = ncol(dat),
  x = as.matrix(dat)
)


# initial values
start_values <- list(
  list("tau"=c(.1, .1, .1, .1, .1, .1, .1),
       "lambda"=c(NA, 0, 0, 0, 0,0,0),
       "inv.phi"=1,
       "inv.psi"=c(1, 1, 1, 1, 1,1,1)),
  list("tau"=c(3, 3, 3, 3, 3,3,3),
       "lambda"=c(NA, 3, 3, 3, 3,3,3),
       "inv.phi"=2,
       "inv.psi"=c(2, 2, 2, 2, 2,2,2)),
  list("tau"=c(5, 5, 5, 5, 5,5,5),
       "lambda"=c(NA, 6, 6, 6, 6,6,6),
       "inv.phi"=.5,
       "inv.psi"=c(.5, .5, .5, .5, .5,.5,.5))
)

# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("tau", "lambda","lambda.std", "phi", "psi")

# fit model
fit <- jags(
  model.file=jags.model,
  data=jags.dat,
  inits=start_values,
  parameters.to.save = param_save,
  n.iter=2000,
  n.burnin = 1000,
  n.chains = 3)

print(fit, width=1000)

