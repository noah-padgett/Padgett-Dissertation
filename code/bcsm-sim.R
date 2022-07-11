# ==================================
# Replicate - Klotzke & Fox (2019)
#               Frontiers Paper
# ==================================
# Created by: R. Noah Padgett
# Created on: 2021-05-17
# Last Edit : 2021-05-17
# ==================================
# Purpose:
# To replicate the results of the
#   simulated dataset from sec 5.
#
# Aim is to reproduce Table 7.
# ==================================

# load packages
source("code/load_packages.R")

# generate data
# sample sizes
nj <- 24 # number of RA+RT
Np <- 200 # number of persons
# matrices
I <- diag(nrow=nj, ncol=nj)
J <- matrix(1, ncol=nj, nrow=nj)
# error variances
sigma <- c(rep(1, nj/2),runif(nj/2, 0.5, 1.5))

# (co)variance parameters
delta <- 0.5
phi <- 0.5
tau <- 0.5
nu <- c(0, -0.05, -0.1, 0.4, 0.2, 0.3)
v <- c(delta, tau, phi, nu)
# u matrix
u <- matrix(
  c(rep(1, 12), rep(0,12),
    rep(0, 12), rep(1,12),
    rep(1, 24),
    rep(c(rep(1, 2), rep(0, 10)),2),
    rep(c(rep(0, 2),rep(1, 2), rep(0, 8)),2),
    rep(c(rep(0, 4),rep(1, 2), rep(0, 6)),2),
    rep(c(rep(0, 6),rep(1, 2), rep(0, 4)),2),
    rep(c(rep(0, 8),rep(1, 2), rep(0, 2)),2),
    rep(c(rep(0, 10), rep(1, 2)),2)),
  ncol=9, nrow=24
)

# expected value no random component
omega <- diag(sigma, nrow=nj, ncol=nj)
i <- 1
for(i in 1:ncol(u)){
  omega<- omega + u[,i,drop=F]%*%t(u[,i,drop=F])*v[i]
}

omega2 <- omega
x <- seq(-5, 5, 0.01)
y <- numeric(length(x))
for(i in 1:length(x)){
  y[i] <- is.positive.definite(omega2 + x[i]*(u[,1]%*%t(u[,2])+u[,2]%*%t(u[,1])))
}
plot(x, y)
dat <- data.frame(x, y)
dat %>%
  filter(y == 1) %>%
  summarise(
    min = min(x), max=max(x)
  )
# data matrix
Y <- rmvnorm(Np, sigma=omega)

# Set up data for JAGS
diagSigma <- matrix(0, ncol=nj, nrow=nj)
diag(diagSigma) <- NA
k <- 9
mydata <- list(
  nj = nj,
  m = Np,
  u=u[,1:k],
  Nt = k,
  Y = Y,
  diagSigma = diagSigma,
  ZeroVec = rep(0, nj)
)

bcsm.model <- function(){

  for(j in 1:m){
    # data model
    Y[j,] ~ dmnorm.vcov(ZeroVec[], omega[,])
  }

  # error variance
  for(i in 1:12){
    # error variance of RT
    sigma[i] ~ dunif(0,3)
    # create diagonal matrix
    diagSigma[i,i] <- 1
    diagSigma[i+12, i+12] <- sigma[i]
  }
  SigmaInv[1, 1:nj, 1:nj] <- inverse(diagSigma[1:nj, 1:nj] )

  # =========================
  # loop over layers of u matrix
  # check value of theta for approx 0 value
  for(t in 1:Nt){
    LAM[t] <- ifelse(abs(theta[t]) < 1e-5, 1e5, 1/theta[t])
    d[t] <- LAM[t] + t(u[,t])%*%SigmaInv[t,1:nj, 1:nj]%*%u[,t]
    # update layer (t) inverse
    A[1:nj, t] <- SigmaInv[t,1:nj, 1:nj]%*%u[1:nj,t]
    # vector multiplication
    for(i in 1:nj){
      for(j in 1:nj){
        B[t,i,j] <- A[i,t]*u[j,t]
      }
    }
    C[t,1:nj,1:nj] <- B[t, , ]%*%SigmaInv[t, , ]
    SigmaInv[t+1, 1:nj, 1:nj] <- SigmaInv[t, 1:nj, 1:nj] - (C[t,,])/d[t]
    # update truncation parameter
    tr[t] <- - 1/(t(u[,t])%*%SigmaInv[t, 1:nj, 1:nj]%*%u[,t])
  }
  # sample theta[t]
  for(t in 1:Nt){
    theta[t] ~ dnorm(0, 0.5);T(tr[t],)
  }
  # compute u%*%t(u)

  for(t in 1:Nt){
    for(i in 1:nj){
      for(j in 1:nj){
        up[t,i,j] <- u[i,t]*u[j,t]
      }
    }
  }
  # add layers of covariance matrix
  for(t in 1:Nt){
    ut[t,1:nj, 1:nj] <- up[t,,]*theta[t]
  }
  for(i in 1:nj){
    for(j in 1:nj){
      utsum[i,j] <- sum(ut[1:Nt,i,j])
    }
  }
  # method interaction
  for(i in 1:nj){
    for(j in 1:nj){
      uprod[i,j] <- u[i,1]*u[j,2] + u[i,2]*u[j,1]
    }
  }
  gamma ~
  # update full covariance matrix
  omega[1:nj , 1:nj] <- diagSigma + utsum

}


# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("sigma", "theta")

# initial values for theta (helps with starting chains)
#   = ensures PD of initial value omega
start.values <- list(
  list("theta"=rep(0.01, k)),
  list("theta"=rep(0.01, k)),
  list("theta"=rep(0.01, k)),
  list("theta"=rep(0.01, k))
)
# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  inits = start.values,
  parameters.to.save = param_save,
  n.iter=2000,n.chains = 4)

print(fit0)
fit <- fit0

jags.mcmc <- as.mcmc(fit)
#jags.mcmc <- as.mcmc(fit$BUGSoutput$sims.matrix)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,regex_pars = "theta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma", prob=0.8)

mcmc_acf(plot.data, regex_pars = "theta")
mcmc_acf(plot.data, regex_pars = "sigma")


mcmc_trace(plot.data, regex_pars = "theta")
mcmc_trace(plot.data, regex_pars = "sigma")

gelman.plot(jags.mcmc)


summary(jags.mcmc)



diagSigma <- matrix(0, ncol=nj, nrow=nj)

k <- 9
nj = nj
m = Np
u=u[,1:k]
Nt = k
Y = Y
theta=rep(0.01, k)
gamma <- -0.1

sigma <- numeric(12)
for(i in 1:12){
  # error variance of RT
  sigma[i] <- runif(1, .5,1.5)
  # create diagonal matrix
  diagSigma[i,i] <- 1
  diagSigma[i+12, i+12] <- sigma[i]
}
SigmaInv <- array(0, dim=c(Nt+2, nj, nj))
tr <- numeric(Nt+2)
SigmaInv[1, , ] <- solve(diagSigma)
t <- 1
for(t in 1:Nt){
  LAM <- ifelse(abs(theta[t]) < 1e-5, 1e5, 1/theta[t])
  d <- LAM + t(u[,t])%*%SigmaInv[t,,]%*%u[,t]
  # update layer (t) inverse
  A <- SigmaInv[t, ,]%*%u[,t]%*%t(u[,t])%*%SigmaInv[t, , ]
  SigmaInv[t+1, , ] <- SigmaInv[t, , ] - (A)/c(d)
  # update truncation parameter
  tr[i] <- - 1/(t(u[,t])%*%SigmaInv[t, , ]%*%u[,t])
}

for(t in 1:Nt){
  x <- -999
  while(x < tr[t]){
    theta[t] = rnorm(1, 0, 2)
    x <- theta[t]
  }
}

utsum <- 0
for(t in 1:Nt){
  utsum <- utsum + u[,t]%*%t(u[,t])*theta[t]
}

# method interaction
d <- LAM + t(u[,1])%*%SigmaInv[Nt+1,,]%*%u[,2]
# update layer (t) inverse
A <- SigmaInv[Nt+1,,]%*%u[,t]%*%t(u[,t])%*%SigmaInv[Nt+1,,]
SigmaInv[Nt+2,,] <- SigmaInv[Nt+1,,] - (A)/c(d)
# update truncation parameter
b1 <- -1/(t(u[,1])%*%SigmaInv[Nt+1,,]%*%u[,2])
b2 <- -1/(t(u[,2])%*%SigmaInv[Nt+1,,]%*%u[,1])

x <- seq(-5, 5, 0.01)
y <- numeric(length(x))
S <- solve(SigmaInv[Nt+1,,])
for(i in 1:length(x)){
  y[i] <- is.positive.definite(S + x[i]*(u[,1]%*%t(u[,2])+u[,2]%*%t(u[,1])))
}
dat <- data.frame(x, y)
dat %>%
  filter(y == 1) %>%
  summarise(
    min = min(x), max=max(x)
  )
