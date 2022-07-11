# Estimate Full Model
# packages needed for this specific file
library(mvtnorm)
library(R2jags)

# Simulate a dataset
# turn into function
simulate_data_misclass <- function(paravec, tau=NULL){
  # NOTE: tau is a matrix[J, C-1] of item threshold parameters that possibly vary over items
  # useful functions
  invlogit <- function(x) {exp(x)/(1+exp(x))}
  logit <- function(x){log(x/(1-x))}
  # Generating Data
  N <- paravec[1] # number of respondents
  J <- paravec[2] # number of items
  C <- paravec[3] # number of response categories
  # ========================= #
  # latent person parameters
  etaCor <- paravec[4] # correlation between ability and speediness
  etasd <- paravec[5:6]
  eta <- mvtnorm::rmvnorm(
    N, mean = c(0, 0),
    sigma = matrix(c(etasd[1], etasd[2]*etaCor,
                     etasd[2]*etaCor, etasd[2]**2),
                   ncol = 2))
  eta0 <- matrix(eta[,1],nrow=1) # ability
  eta1 <- matrix(eta[,2],nrow=1) # log speediness
  # ========================= #
  # item parameters
  # item factor loadings
  lambda <- matrix(rep(paravec[7], J), ncol=1)
  # item latent residual variances
  theta <- c(1 - lambda**2)
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
  ystar <- apply(ystar, 2, FUN = function(x){mvtnorm::rmvnorm(1, x, diag(theta, ncol=J, nrow=J))})
  # response time parameters (copied from Molenaar et al. 2021)
  nu <- matrix(rep(paravec[8], J), ncol=1)
  sigma.ei <- matrix(rep(paravec[9], J), ncol=1)
  rho1 <- paravec[10]
  #rho2 <- 0
  #delta <- 0

  mulogt <- logt <- matrix(nrow=N, ncol=J)
  i<-j <- 1
  for(i in 1:N){
    for(j in 1:J){
      # obtain expected log response time
      mulogt[i,j] <- nu[j, 1] - eta1[1,i] - rho1*abs( eta0[1,i] - sum(tau[j,])/length(tau[j,]) )
      # sample observed log response time
      # logRT ~ N(mulogt, sigma.ie)
      logt[i,j] <- rnorm(1, mulogt[i,j], sqrt(sigma.ei[j,1]))
    }
  }

  # construct missclassification
  # based on latent response time (nu - eta1)
  misclass.time.trans <- function(lrt, c, b, K, diagonal = FALSE){
    if(c == b){
      g <- 1/(1 + exp(-lrt))
      if(diagonal == TRUE){
        g <- 1
      }
    }
    if(c != b){
      g <- (1/(K-1))*(1-1/(1 + exp(-lrt)))
      if(diagonal == TRUE){
        g <- 0
      }
    }

    g

  }

  gamma <- array(dim=c(N,J,C,C))

  for(i in 1:N){for(j in 1:J){for(b in 1:C){for(c in 1:C){
    gamma[i,j,b,c] <- misclass.time.trans(nu[j, 1] - eta1[1, i], b, c, C)
  }}}}# end loops


  pi <- pi.gte <- omega <- array(0,dim=c(N, J, C))
  Y <- matrix(nrow=N, ncol=J)
  i <- j <- c <- 1
  for(i in 1:N){
    for(j in 1:J){

      # GRM model
      for(c in 2:C){
        # P(greater than or equal to category c > 1)
        pi.gte[i,j,c] <- invlogit(ystar[j,i]-tau[j,(c-1)])
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

      # observed category prob (Pr(y=c))
      for(c in 1:C){
        for(ct in 1:C){
          # sum over ct
          omega[i,j,c] = omega[i,j,c] + gamma[i,j,ct,c]*pi[i,j,ct]
        }
      }
      Y[i,j] <- sample(x=1:C, size=1, prob=omega[i,j,])
    }
  }
  # true_values <- list(eta0, eta1, lambda, nu, sigma.ei, tau, mulogt, ystar, theta, gamma, omega)
  # names(true_values) <- c("eta", "")
  sim_data <- list(Y, logt)
  names(sim_data) <- c("y", "logt")
  return(sim_data)

}

# get data
paravec <- c(
  N = 500, J = 5, C = 5,
  etaCor = .23, etasd1 = 1, etasd2 = sqrt(0.1),
  lambda=0.9, nu=1.5, sigma.ei=0.25, rho1=0.1)
sTau <- matrix(
  c(-2.5, -0.5, 1, 1.5,
    -2, -0.5, 1.5, 2.25,
    -2.25, -0.5, 1.5, 2.1,
    -1.75, -0.5, 2, 2.5,
    -1.75, -0.5, 1.75, 2.5),
  ncol=paravec[3]-1, nrow=paravec[2], byrow=T
)
sim.data <- simulate_data_misclass(paravec, tau=sTau)

mydata <- list(
  y = sim.data$y,
  N = nrow(sim.data$y),
  nit=ncol(sim.data$y)
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # data model
      y[p,i] ~ dcat(pi[p,i, ])

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

    }
  }
  ### Priors
  # person parameters
  for(p in 1:N){
    eta[p] ~ dnorm(0, 1) # latent ability
  }
  lambda[1] = 1
  theta[1] = 1 + pow(lambda[1],2)
  tau[1, 1] = 0
  tau[1, 2] ~ dnorm(0, 0.1);T(tau[1, 1],)
  tau[1, 3] ~ dnorm(0, 0.1);T(tau[1, 2],)
  tau[1, 4] ~ dnorm(0, 0.1);T(tau[1, 3],)
  for(i in 2:nit){
    # Thresholds
    tau[i, 1] = 0
    tau[i, 2] ~ dnorm(0, 0.1);T(tau[i, 1],)
    tau[i, 3] ~ dnorm(0, 0.1);T(tau[i, 2],)
    tau[i, 4] ~ dnorm(0, 0.1);T(tau[i, 3],)
    # loadings
    lambda[i] ~ dnorm(0, .44);T(0,)
    # LRV total variance
    # total variance = residual variance + fact. Var.
    theta[i] = 1 + pow(lambda[i],2)
  }
  # omega
  l2 = (lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5])
  reli.omega = (pow(l2,2))/(pow(l2,2)+ nit)
}

jags.params <- c(
  # ability measurement model
  "tau", "lambda", "theta", "pi",
  # reliability
  "reli.omega"
)
jags.m <- rjags::jags.model(
  file="code/kodiak-files/jags_model.txt",
  data = mydata
)

jags.m$recompile()
jags.m.samples <- jags.samples(
  jags.m, jags.params,
  n.iter=1000)
jags.m.samples

# fit model
a <- Sys.time()
fit <- jags.parallel(
  model.file=jags.model
  , data=mydata
  , parameters.to.save = jags.params
  , n.iter=5000
  , n.chains = 2
)
b <- Sys.time()
b-a
fit.summary <- fit$BUGSoutput$summary

