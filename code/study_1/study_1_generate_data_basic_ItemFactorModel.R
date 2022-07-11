# ===================================== #
# study_1_generate_data.R
# ===================================== #
# Padgett - Dissertation
# Created
#   on: 2022-01-06
#   by: R. Noah Padgett
# Last Editted
#   on: 2022-01-07
#   by: R. Noah Padgett
# ===================================== #
# Purpose: Generate data for study 1
# ===================================== #
# libraries
library(R2jags)
library(mvtnorm)

set.seed(1)
N <- 100
N_cat <- 3
N_items <- 5

# data parameters
paravec <- c(
	  N = N
	, J = N_items
	, C = N_cat
	, etaCor = .23
	, etasd1 = 1
	, etasd2 = sqrt(0.1)
  , lambda=0.9
	, nu=1.5
	, sigma.ei=0.25
	, rho1=0.1
)
# thresholds
sim_tau <- matrix(ncol=N_cat-1, nrow=N_items)
for(c in 1:(N_cat-1)){
  if(c == 1){
    sim_tau[,1] <- runif(N_items, -1, -0.33)
  }
  if(c > 1){
    sim_tau[,c] <- sim_tau[,c-1] + runif(N_items, 0.33, 1.33)
  }
}


simulate_data_ifa <- function(paravec, tau=NULL, scaling = "marginal", is.probit=F){
  invlogit <- function(x) {exp(x)/(1+exp(x))}
  logit <- function(x){log(x/(1-x))}

  invprobit <- function(x) {pnorm(x)}
  probit <- function(x){qnorm(x)}

  # Generating Data
  N <- paravec[1] # number of respondents
  J <- paravec[2] # number of items
  C <- paravec[3] # number of response categories
  # ========================= #
  # latent person parameters
  etaCor <- paravec[4] # correlation between ability and speediness
  etasd <- paravec[5:6]
  eta <- rnorm(N, 0, 1)
  eta0 <- matrix(eta, nrow=1) # ability
  # ========================= #
  # item parameters
  # item factor loadings
  lambda <- matrix(rep(paravec[7], J), ncol=1)
  # latent response scaling
  # fix residual of latent response to 1
  #   i.e. probit regression
  if(scaling == "conditional"){
    theta = matrix(rep(1, J), ncol=1)
  }
  # fix total variance of latent response to 1
  if(scaling == "marginal"){
    theta = 1 - lambda**2*etasd[1]
  }


  # item thresholds
  if(is.null(tau)){
    tau <- matrix(ncol=C-1, nrow=J)
    for(c in 1:(C-1)){
      if(c == 1){
        tau[,1] <- runif(J, -1, -0.33)
      }
      if(c > 1){
        tau[,c] <- tau[,c-1] + runif(J, 0.25, 1)
      }
    }
  }

  # latent item response
  ystar <- lambda%*%eta0
  for(i in 1:J){
    for(p in 1:N){
      ystar[i,p] <- rnorm(1, ystar[i,p], theta[i,1])
    }
  }

  pi <- pi.gte <- omega <- array(0,dim=c(N, J, C))
  Y <- matrix(nrow=N, ncol=J)
  i <- j <- c <- 1
  for(i in 1:N){
    for(j in 1:J){

      # GRM model
      for(c in 2:C){
        # P(greater than or equal to category c > 1)
        if(is.probit == F){
          pi.gte[i,j,c] <- invlogit(ystar[j,i]-tau[j,(c-1)])
        }
        if(is.probit == T){
          pi.gte[i,j,c] <- invprobit(ystar[j,i]-tau[j,(c-1)])
        }
      }
      # P(greater than or equal to category 1)
      pi.gte[i,j,1] <- 1
      # equal to prob.
      for(c in 1:(C-1)){
        # P(greater equal to category c < C)
        pi[i,j,c] <- pi.gte[i,j,c]-pi.gte[i,j,c+1]
      }
      # P(greater equal to category C)
      pi[i,j,C] <- pi.gte[i,j,C]


      Y[i,j] <- sample(x=1:C, size=1, prob=pi[i,j,])
      # rescale to 0/1 if dichotomous items
      if(C == 2){
        Y[i,j] = Y[i,j]-1
      }
    }
  }
  # true_values <- list(eta0, eta1, lambda, nu, sigma.ei, tau, mulogt, ystar, theta, gamma, omega)
  sim_data <- list(Y)
  names(sim_data) <- c("y")

  return(sim_data)

}

# Use parameters to simulate data
sim.data.marg <- simulate_data_ifa(paravec, tau=sim_tau, is.probit=T)
sim.data.cond <- simulate_data_ifa(paravec, tau=sim_tau, is.probit=T, scaling="conditional")


# Save parameters
jags.params <- c(
  "tau", "lambda", "theta", "reli.omega"
)

# marginal data
mydata.marg <- list(
  y = sim.data.marg$y,
  N = 100,
  nit=5
)
# conditional data
mydata.cond <- list(
  y = sim.data.cond$y,
  N = 100,
  nit=5
)

# parameterization 1
D_marg.M_marg <-  R2jags::jags(
  model=paste0(getwd(),"/code/study_1/model_1_marg.txt")
  , parameters.to.save = jags.params
  , data=mydata.marg
  , n.chains = 4
  , n.burnin = 1000
  , n.iter = 2000
)
print(D_marg.M_marg)
# parameterization 2
D_marg.M_cond <-  R2jags::jags(
  model=paste0(getwd(),"/code/study_1/model_1_cond.txt")
  , parameters.to.save = jags.params
  , data=mydata.marg
  , n.chains = 4
  , n.burnin = 1000
  , n.iter = 2000
)
print(D_marg.M_cond)
# parameterization 3
D_cond.M_marg <-  R2jags::jags(
  model=paste0(getwd(),"/code/study_1/model_1_marg.txt")
  , parameters.to.save = jags.params
  , data=mydata.cond
  , n.chains = 4
  , n.burnin = 1000
  , n.iter = 2000
)
print(D_cond.M_marg)
# parameterization 4
D_cond.M_cond <-  R2jags::jags(
  model=paste0(getwd(),"/code/study_1/model_1_cond.txt")
  , parameters.to.save = jags.params
  , data=mydata.cond
  , n.chains = 4
  , n.burnin = 1000
  , n.iter = 2000
)
print(D_cond.M_cond)


