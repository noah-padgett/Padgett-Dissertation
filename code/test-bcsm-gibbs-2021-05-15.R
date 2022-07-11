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
library(mvtnorm)
library(LaplacesDemon)
library(ggplot2)
library(patchwork)


# truncated shifted inverse-gamma distribution
dIG <- function(x, alpha0, beta0, psit=0, trt=0){
  # check if x + psi > 0
  xshift <- ifelse(x + psit > 0, x+psit, 1e-5)
  out <- ((beta0^alpha0)/gamma(alpha0))*(xshift)^(-alpha0 -1)*exp(-beta0/(xshift))
  for(i in 1:length(out)){
    out[i] <- ifelse(out[i] < trt, 0, out[i])
  }
  out
}
pIG <- function(x0, alpha0, beta0, psit=0, trt=0, ll=0,...){
  integrate(dIG, ll, x0,...)$value
}

# this works!!!!
sliceSample = function (n, f, x.interval = c(0, 1), root.accuracy = 0.01, root.wiggle=1e-5) {
  # n <- 1
  # f <- target
  # x.interval <- c(0, upperBound)
  # root.accuracy = upperBound/1000
  # root.wiggle <- 1e-5
  # n is the number of points wanted
  # f is the target distribution
  # x.interval is the A,B range of x values possible.
  pts = vector("numeric", n) # This vector will hold our points.
  x = runif(1, x.interval[1], x.interval[2]) # Take a random starting x value.
  i = 1
  for (i in 1:n) {
    pts[i] = x
    y = runif(1, 0, f(x)) # Take a random y value
    # Imagine a horizontal line across the distribution.
    # Find intersections across that line.
    fshift = function (x) { f(x) - y }
    roots = c()
    j = x.interval[1]
    for (j in seq(x.interval[1]+root.wiggle, x.interval[2] - root.accuracy, by = root.accuracy)) {
      if ((fshift(j) < 0) != (fshift(j + root.accuracy) < 0)) {
        # Signs don't match, so we have a root.
        root = uniroot(fshift, c(j, j + root.accuracy))$root
        roots = c(roots, root)
      }
    }
    # Include the endpoints of the interval.
    roots = c(x.interval[1], roots, x.interval[2])
    # Divide intersections into line segments.
    segments = matrix(ncol = 2)
    for (j in 1:(length(roots) - 1)) {
      midpoint = (roots[j + 1] + roots[j]) / 2.0
      if (f(midpoint) > y) {
        # Since this line segment is under the curve, add it to segments.
        segments = rbind(segments, c(roots[j], roots[j + 1]))
      }
    }
    # Uniformly sample next x from segments
    # Assign each segment a probability, then unif based on those probabilities.
    # This is a bit of a hack because segments first row is NA, NA
    # Yet, subsetting it runs the risk of reducing the matrix to a vector in special case.
    total = sum(sapply(2:nrow(segments), function (i) {
      segments[i, 2] - segments[i, 1]
    }))
    probs = sapply(2:nrow(segments), function (i) {
      (segments[i, 2] - segments[i, 1]) / total
    })
    # Assign probabilities to each line segment based on how long it is
    # Select a line segment by index (named seg)
    p = runif(1, 0, 1)
    selectSegment = function (x, i) {
      if (p < x) return(i)
      else return(selectSegment(x + probs[i + 1], i + 1))
    }
    seg = selectSegment(probs[1], 1)
    # Uniformly sample new x value
    x = runif(1, segments[seg + 1, 1], segments[seg + 1, 2])
  }
  return(pts)
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
lambda <- matrix(rnorm(nItem, 0, 1), ncol=1, nrow=nItem)
Sigma0 <- diag(1, ncol=nItem, nrow=nItem)
# covariance parameters
Delta <- c(0)
# Design Matrix
u <- matrix(
  c(1, 1, 1, 1, 1),
  nrow=nItem, ncol=1
)
Nt <- ncol(u)
# covariance structure
Sigma <- Sigma0 + Delta*(u%*%t(u))

# sample data
Y <- matrix(0, ncol = nItem, nrow=Np)
for(p in 1:Np){
  Y[p, 1:nItem] <- t(lambda) + mvtnorm::rmvnorm(1, sigma = Sigma)
}
# center to remove mean structure (for now)
Y <- apply(Y, 2, function(x) { x - mean(x) })


# parameters to monitor
Niter <- 1000
sigma <- matrix(nrow=Niter, ncol=nItem)
SigmaFull <- array(dim=c(Niter, Nt+1, nItem, nItem))
theta <- matrix(nrow=Niter, ncol=1)
psi <- matrix(nrow=Niter, ncol=1)
trt <- matrix(nrow=Niter, ncol=1)
# initial values
sigma[1,] <- runif(nItem, 0.5, 1.5) # item unique variances
SigmaFull[1,1,,] <- diag(sigma[1,], ncol=nItem, nrow=nItem) # full covariance matrix at different layers
theta[1,] <- 0 # covariance parameter
psi[1,] <- 0
trt[1,] <- 0

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
  #   Compute SSW(k)
  Y. <- colMeans(Y)
  Y.mat <- matrix(
    rep(Y., nrow(Y)),
    nrow=nrow(Y),
    ncol=length(Y.),
    byrow = T
  )
  SSW <- t(Y - Y.mat)%*%(Y - Y.mat)

  # ====================
  #   Sample sigma[i]
  # make it constant across items
  SSWi <- mean(diag(SSW))
  Q <- theta[iter-1,t]*u[i,t]
  sigma[iter, i] <- rTSIG(1, alpha0 + nItem/2, beta0 + SSBt/2, Q, 0)
  # target <- function(x){ dIG(x, alpha0 + nItem/2, beta0 + SSWi/2, Q, 0) }
  # upperBound <- 1/qgamma(.001, alpha0 + nItem/2, beta0 + SSWi/2)
  # sigma[iter, i] = sliceSample(n = 1, target, x.interval=c(0, upperBound), root.accuracy = upperBound/1000)

  # for(i in 1:nItem){
  #   SSWi <- SSW[i,i]
  #   Q <- theta[iter-1,t]*u[i,t]
  #   target <- function(x){ dIG(x, alpha0 + nItem/2, beta0 + SSWi/2, Q, 0) }
  #   upperBound <- 1/qgamma(.001, alpha0 + nItem/2, beta0 + SSWi/2)
  #   sigma[iter, i] = sliceSample(n = 1, target, x.interval=c(0, upperBound), root.accuracy = upperBound/1000)
  # }

  # Computer SSB eq (8)
  Yut <- Y[,which(u[,t] != 0)]
  Y.t <- matrix(rep(mean(Yut), nrow(Yut)), ncol=1)
  Y.it <- matrix(rowMeans(Yut), ncol=1)
  SSBt <- c(t(Y.it - Y.t)%*%(Y.it - Y.t))
  # =====================
  # # truncation parameter
  # SigmaFull[iter,1, , ] <- diag(sigma[iter, ], ncol=nItem, nrow=nItem)
  # SigmaFull[iter, t+1, , ] <- diag(sigma[iter, ], ncol=nItem, nrow=nItem) + theta[iter-1, t] * (u[,t, drop=F]%*%t(u[,t, drop=F]))
  # # get iterative inverses
  # SFt <- SigmaFull[iter,t, , ]
  # SFtp1 <- SigmaFull[iter,t+1, , ]
  # # inverse of full matrix (layer t=0)
  # SF0inv <- solve(SigmaFull[iter, 1, , ])
  # if(t == 1){
  #   SFtinv <- SF0inv
  # } else {
  #   # compute denominator of eq (5)
  #   invtheta <- ifelse(theta[iter-1,t] == 0, 0, 1/theta[iter-1,t])
  #   devisor <- c(invtheta + t(u[,t, drop=F])%*%SF0inv%*%u[,t, drop=F])
  #   # compute Inverse at layer t-1
  #   SFtinv <- SF0inv - (SF0inv%*%u[,t, drop=F]%*%t(u[,t, drop=F])%*%SF0inv)/( devisor )
  # }

  # compute lower bound of theta[t]
  #trt[iter, t] <- -1/(t(u[,t, drop=F])%*%SFtinv%*%u[,t, drop=F])
  trt[iter, t] <- - sigma[iter,i]/nItem
  if(trt[iter, t] > 0) trt[iter, t] <- 0
  # =====================
  # shift parameter
  psi[iter, t] <- SSBt - theta[iter-1, t]
  # =====================
  #   Sample theta
  # sample the parameter from posterior
  theta[iter, t] <- rTSIG(1, alpha0 + N/2, beta0 + SSBt/2, psi[iter, t], trt[iter, t])
  #
  if(iter %% 100 == 0) cat(".")
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


a <- rTSIG(1000, alpha0 + N/2, beta0 + SSBt/2, psi[iter,t], trt[iter, t])
b <- rinvgamma(1000, alpha0 + N/2, beta0 + SSBt/2)
plot.dat <- data.frame(A=a, B=b)
ggplot(plot.dat)+
  geom_density(aes(x=A), color="black")+
  geom_density(aes(x=B), color="blue")+
  theme_classic()

mcmc.draws <- cbind(theta, sigma[,1], psi, trt)
colnames(mcmc.draws) <- c("theta[1]", paste0("sigma[",1,"]"), "pst", "trt")

library(coda)
library(bayesplot)
mcmc.draws <- as.mcmc(mcmc.draws)
summary(mcmc.draws)

mcmc_areas(mcmc.draws)
mcmc_acf(mcmc.draws)

traceplot(mcmc.draws)

