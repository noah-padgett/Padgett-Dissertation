# ================================== #
# Attempting Gibbs Sampler for BCSM
#
# create by: R. Noah Padgett
# create on: 2021-05-11
# last edit: 2021-05-15
# ================================== #
# Reference:
# Klotzke & Fox (2019a). BCSM. Psychometrika.
# Klotzke & Fox (2019b). BCSM of RA and RT. Frontiers.
# ================================== #

# load packages
source("code/load_packages.R")


# truncated shifted inverse-gamma distribution
dIG <- function(x, alpha0, beta0, psit=0, trt=0){
  # check if x + psi > 0
  xshift <- ifelse(x + psit > 0, x+psit, 0)
  out <- ((beta0^alpha0)/gamma(alpha0))*(xshift)^(-alpha0 -1)*exp(-beta0/(xshift))
  for(i in 1:length(out)){
    out[i] <- ifelse(out[i] < trt, 0, out[i])
  }
  out
}
pIG <- function(x0, alpha0, beta0, psit=0, trt=0, ll=0,...){
  integrate(dIG, ll, x0,...)$value
}
rTSIG <- function(n, alpha0, beta0, psi0, trt0, max.iter=1000){
  x <- numeric(n)
  for(j in 1:n){
    i <- 1
    x[j] <- trt0
    while(x[j] <= trt0){
      x[j] <- 1/rgamma(1, alpha0, beta0)
      x[j] <- x[j] - psi0
      i <- i + 1
      if(i > max.iter) break
    }
  }
  x
}

# updating functions
SigmaTinv.update <- function(Atinv, lambda, nu){
  LAM <- ifelse(abs(lambda) < 1e-5, 1e5, 1/lambda)
  d <- c(LAM + t(nu)%*%Atinv%*%u)
  Atinv - (Atinv%*%nu%*%t(nu)%*%Atinv)/d
}

# Simulate date for test
Np <- 50 # num persons
nItem <- 5 # number items
N <- Np*nItem
# person matrices
Ip <- diag(1, ncol=Np, nrow=Np)
Jp <- matrix(1, ncol=Np, nrow=Np)
# item matrices
Ii <- diag(1, ncol=nItem, nrow=nItem)
Ji <- matrix(1, ncol=nItem, nrow=nItem)
# other helpful matrices
Jvec <- matrix(1, ncol=1, nrow=nItem)
# item parameters
mu <- matrix(rnorm(nItem, 0, 1), ncol=1, nrow=nItem)
sig0 <- 1
# covariance parameters
theta0 <- -sig0/nItem + 0.05
# Design Matrix
u <- matrix(
  c(1, 1, 1, 1, 1),
  nrow=nItem, ncol=1
)
Nt <- ncol(u)
# covariance structure
Sigma <- Ii*sig0 + theta0*(u%*%t(u))

# sample data
Y <- matrix(0, ncol = nItem, nrow=Np)
for(p in 1:Np){
  Y[p, 1:nItem] <- t(mu) + mvtnorm::rmvnorm(1, sigma = Sigma)
}
# center to remove mean structure (for now)
Y <- apply(Y, 2, function(x) { x - mean(x) })


# parameters to monitor
Niter <- 5000
sigma <- matrix(nrow=Niter, ncol=nItem)
SigmaFull <- array(dim=c(Niter, Nt+1, nItem, nItem))
theta <- matrix(nrow=Niter, ncol=1)
psi <- matrix(nrow=Niter, ncol=1)
trt <- matrix(nrow=Niter, ncol=1)
# initial values
sigma[1,] <- rep(1, nItem) # residual varaince
SigmaFull[1,1,,] <- diag(sigma[1,], ncol=nItem, nrow=nItem) # full covariance matrix at different layers
theta[1,] <- 0 # covariance parameter
psi[1,] <- 0 # shift parameter
trt[1,] <- 0 # truncation parameter

# priors
alpha0 <- 0.1
beta0 <- 0.1
# iterators
i <- t <- 1
iter <- 2
# loop
for(iter in 2:Niter){
  i <- t <- 1
  # ====================
  #   Compute SSB[t]
  SSB <- numeric(Nt)
  Yp.ut <- rowMeans(Y[, which(u[,t] != 0)])
  Y..ut <- mean(Y[, which(u[,t] != 0)])
  SSB[t] <- sum((Yp.ut - Y..ut)**2)

  # ====================
  #   Compute SSW[k]
  SSW <- numeric(nItem)
  Yp. <- colMeans(Y)
  for(i in 1:nItem){
      SSW[i] <- sum((Y[,i] - Yp.[i])**2)
  }

  # ====================
  #   Sample sigma[i]
  for(i in 1:nItem){
    Q <- c(theta[iter-1,t]*u[i,t])
    sigma[iter, i] <- rTSIG(1, alpha0 + nItem/2, beta0 + SSB[t]/2, Q, 0)
  }

  # ====================
  #   Compute truncation
  SigmaTinv <- array(0,c(Nt+1, nItem, nItem))
  # base layer Sigma0inv
  SigmaTinv[1, , ] <- solve(diag(sigma[iter,], ncol=nItem, nrow=nItem))
  # next layer
  SigmaTinv[2, , ] <- SigmaTinv.update(SigmaTinv[1, , , drop=T], theta[iter-1,1], u[,1, drop=F])

  # =====================
  # # truncation parameter
  trt[iter, t] <- - 1/(t(u[,t,drop=F])%*%SigmaTinv[1, , , drop=T]%*%u[,t,drop=F])
  if(trt[iter, t] > 0) trt[iter, t] <- 0
  # =====================
  # shift parameter
  # shift is the negative of the truncation parameter
  # Fox et al. (2021, p. 226)
  psi[iter, t] <- - trt[iter, t]
  # =====================
  #   Sample theta
  # sample the parameter from posterior
  theta[iter, t] <- rTSIG(1, alpha0 + nItem/2, beta0 + SSB[t]/2, psi[iter, t], trt[iter, t])
  #

  if(iter == nBurn) iter <- 1
  if(iter %% 100 == 0) cat(".")

}




# structure draws
mcmc.draws <- cbind(theta, sigma, psi, trt)
mcmc.draws <- mcmc.draws[-c(1:2500),]
colnames(mcmc.draws) <- c("theta[1]", paste0("sigma[",1:nItem,"]"), "pst", "trt")

mcmc.draws <- as.mcmc(mcmc.draws)
summary(mcmc.draws)

# remove burn-in
mcmc.draws <- mcmc.draws[-c(1:999, seq(1000,5000,20)),]
mcmc.draws <- as.mcmc(mcmc.draws)
summary(mcmc.draws)

mcmc_areas(mcmc.draws)
mcmc_acf(mcmc.draws)

traceplot(mcmc.draws)

