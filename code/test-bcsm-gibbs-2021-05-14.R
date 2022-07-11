# ================================== #
# Attempting Gibbs Sampler for BCSM
#
# create by: R. Noah Padgett
# create on: 2021-05-11
# last edit: 2021-05-11
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

p <- 100 # num persons
nItem <- 10 # number items
# person matrices
Ip <- diag(1, ncol=p, nrow=p)
Jp <- matrix(1, ncol=p, nrow=p)
# item matrices
Ii <- diag(1, ncol=nItem, nrow=nItem)
Ji <- matrix(1, ncol=nItem, nrow=nItem)
# other helpful matrices
Jvec <- matrix(1, ncol=1, nrow=nItem)
# item parameters
lambda <- matrix(rnorm(nItem, 0, 1), ncol=1, nrow=nItem)
Sigma0 <- diag(1, ncol=nItem, nrow=nItem)
# covariance parameters
Delta <- c(0, 0.01, 0.05)
Nt <- length(Delta)
delta <- 0.2

# Design Matrix
u <- matrix(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 0, 1, 1, 1, 1, 0, 0, 0),
  nrow=nItem, ncol=3
)


# covariance structure
Sigma <- Sigma0 + delta*Ji
for(j in 1:Nt){
  Sigma <- Sigma + Delta[j]*(u[,j, drop=F]%*%t(u[,j, drop=F]))
}

# sample data
Y <- matrix(0, ncol = nItem, nrow=p)
for(i in 1:p){
  Y[i, 1:nItem] <- t(lambda) + mvtnorm::rmvnorm(1, sigma = Sigma)
}

# base layer of (co)variance matrix
deltaMin <- -1/(t(Jvec)%*%solve(Sigma0)%*%Jvec)
deltaMin
deltaUse <- as.numeric(ifelse(delta > deltaMin, delta, deltaMin+0.0001))
A1 <- Sigma0 + deltaUse*Ji

Ad <- A1

psi <- 1
nu <- u[,1, drop=F]

Adp1 <- Ad + psi*(nu%*%t(nu))

psiMin <- -1/(t(nu)%*%solve(Ad)%*%nu)


# truncated shifted inverse-gamma distribution
dIG <- function(x, alpha0, beta0, psit=0, trt=0){
  ((beta0^alpha0)/gamma(alpha0))*(x+psit)^(-alpha0 -1)*exp(-beta0/(x+psit))
}
pIG <- function(x0, alpha0, beta0, psit=0, trt=0, ll=0,...){
  integrate(dIG, ll, x0,...)$value
}


# this works!!!!
sliceSample = function (n, f, x.interval = c(0, 1), root.accuracy = 0.01, root.wiggle=1e-5) {
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


target <- function(x){ dIG(x, 10, 1, 5, 0) }


points = sliceSample(n = 1000, target, x.interval=c(0, 5), root.accuracy = 0.001)


plot.dat1 <- data.frame(samples=points)
plot.dat2 <- data.frame(X=seq(-0.05, 1, 0.001), Y=target(seq(-0.05, 1, 0.001)))

p1 <- ggplot(plot.dat1, aes(x=samples))+
  geom_density()+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
p2 <- ggplot(plot.dat2, aes(x=X, y=Y))+
  geom_line()+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
p1+p2

# computing SSB (8) - sufficient stat for (7)
k <- t <- i <-  1

# identify which items are in vector u[t]
Yut <- Y[,which(u[,t] != 0)]
# item means
Y.t <- mean(Yut)
# person means
Y.it <- rowMeans(Yut)
# sum of squares of deviations of the conditional person level means from the conditional grand mean
# Sum of squares for effect u[t]
SSBt <- sum((Y.it - Y.t)**2)
# OR equivalently
Y.t <- matrix(rep(mean(Yut), nrow(Yut)), ncol=1)
Y.it <- matrix(rowMeans(Yut), ncol=1)
t(Y.it - Y.t)%*%(Y.it - Y.t)


# compute SSW (9)
Y. <- colMeans(Y)
Y.mat <- matrix(
  rep(Y., nrow(Y)),
  nrow=nrow(Y),
  ncol=length(Y.),
  byrow = T
)
SSW <- t(Y - Y.mat)%*%(Y - Y.mat)
SSWk <- SSW[k,k]


# compute psit (shift parameters)
#psit <- SSB[t] - theta[t]





# Simulate date for test
Np <- 50 # num persons
nItem <- 5 # number items
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
Delta <- c(0.1)
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
Niter <- 100
sigma <- matrix(nrow=Niter, ncol=nItem)
Sigma <- array(dim=c(Niter, nItem, nItem))
SigmaFull <- array(dim=c(Niter, Nt+1, nItem, nItem))
Sigma11 <- array(dim=c(Niter, nItem, 1, 1))
Sigma12 <- array(dim=c(Niter, nItem, 1, nItem-1))
Sigma21 <- array(dim=c(Niter, nItem, nItem-1, 1))\
Sigma22 <- array(dim=c(Niter, nItem, nItem-1, nItem-1))
Sigma0inv <- array(dim=c(Niter, nItem,  nItem-1, nItem-1))
Sigmatinv <- array(dim=c(Niter, Nt, nItem, nItem-1, nItem-1))
theta <- matrix(nrow=Niter, ncol=1)
condSigma2 <- matrix(nrow=Niter, ncol=nItem)

# initial values
sigma[1,] <- runif(nItem, 0.5, 1.5) # item unique variances
Sigma[1,,] <- diag(sigma[1,], ncol=nItem, nrow=nItem) # initialize covariance matrix
SigmaFull[1,1,,] <- diag(sigma[1,], ncol=nItem, nrow=nItem) # full covariance matrix at different layers
Sigma11[1,1, ,] <- Sigma[1,1,1]
Sigma12[1,1,,] <- Sigma[1,1,-1]
Sigma21[1,1,,] <- Sigma[1,-1,1]
Sigma22[1,1,,] <- Sigma[1,-1,-1]
Sigma0inv[1,1,,] <- solve(Sigma[1,-1,-1])
Sigmatinv[1,1,,,] <- solve(Sigma[1,-1,-1])
theta[1,] <- 0 # covariance parameter


# priors
alpha0 <- 0.1
beta0 <- 0.1


iter <- i <- t <- 1

# create conditional variances
for(i in 1:nItem){
  # extract partitioned covariance matrix
  # decomposition of Sigma
  S11 <- Sigma11[iter,i, , ] <- Sigma[iter, i, i, drop=F]
  S12 <- Sigma12[iter,i, , ] <- Sigma[iter, i, -i, drop=F]
  S21 <- Sigma21[iter,i, , ] <- Sigma[iter, -i, i, drop=F]
  S22 <- Sigma22[iter,i, , ] <- Sigma[iter, -i, -i, drop=F]
  S11 <- S11[ , , , drop=T]
  S12 <- matrix(S12[,,,drop=T], nrow=1, ncol=nItem-1)
  S21 <- matrix(S21[,,,drop=T], nrow=nItem-1, ncol=1)
  S22 <- matrix(S22[,,,drop=T], nrow=nItem-1, ncol=nItem-1)
  S22inv <- solve(S22)


  condSigma2[iter, i] <- c(S11) - c(S12%*%S22inv%*%S21)

}

# design matrix layer
ut <- u[-i,t, drop=F]
# covariance parameter check for 0
invtheta <- ifelse(theta[t] == 0, 0, 1/theta[t])
devisor <- c(invtheta + t(ut)%*%S0inv%*%ut)
Sigmatinv[iter, t, i, , ] <- S0inv - (S0inv%*%ut%*%t(ut)%*%S0inv)/( devisor )
Stinv <- Sigmatinv[iter, t, i, , ]
# obtain conditional variance
condSigma2[iter, i] <- c(S11) - c(S12%*%Stinv%*%S21)





# sample data # data model - not needed?
Yiter <- matrix(nrow=Np, ncol=nItem)
for(p in 1:Np){
  for(i in 1:nItem){
    Yiter[p, i] <- rnorm(1, 0, condSigma2[iter, i])
  }
}


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
SSWi <- SSW[i,i]
# ====================
#   Sample sigma[i]
Q <- theta[iter,t]*u[i,t]
target <- function(x){ dIG(x, alpha0 + Np/2, beta0 + SSWi/2, Q, 0) }
upperBound <- 1/qgamma(.001, alpha0 + Np/2, beta0 + SSWi/2)
sigma[iter, i] = sliceSample(n = 1, target, x.interval=c(0, upperBound), root.accuracy = upperBound/1000)
# Computer SSB eq (8)
Yut <- Y[,which(u[,t] != 0)]
Y.t <- matrix(rep(mean(Yut), nrow(Yut)), ncol=1)
Y.it <- matrix(rowMeans(Yut), ncol=1)
SSBt <- c(t(Y.it - Y.t)%*%(Y.it - Y.t))
# =====================
# shift parameter
psi <- SSBt - theta[t]
# =====================
# truncation parameter
SigmaFull[iter,1, , ] <- diag(sigma[iter, ], ncol=nItem, nrow=nItem)
SigmaFull[iter, t+1, , ] <- diag(sigma[iter, ], ncol=nItem, nrow=nItem) + theta[iter, t] * (u[,t, drop=F]%*%t(u[,t, drop=F]))
# get iterative inverses
SFt <- SigmaFull[iter,t, , ]
SFtp1 <- SigmaFull[iter,t+1, , ]
# inverse of full matrix
SF0inv <- solve(SigmaFull[iter, 1, , ])
# compute denominator of eq (5)
invtheta <- ifelse(theta[t] == 0, 0, 1/theta[t])
devisor <- c(invtheta + t(u[,t, drop=F])%*%SF0inv%*%u[,t, drop=F])
# compute Inverse at layer t
SFtinv <- SF0inv - (SF0inv%*%u[,t, drop=F]%*%t(u[,t, drop=F])%*%SF0inv)/( devisor )

# compute lower bound of theta[t]
trt <- -1/(t(u[,t, drop=F])%*%SFtinv%*%u[,t, drop=F])
if(trt > 0) trt <- 0

# =====================
#   Sample theta
# set up posterior
target <- function(x){ dIG(x, alpha0 + Np/2, beta0 + SSBt/2, psi, trt) }
# compute an upper bound for slice sampling
upperBound <- 1/qgamma(.001, alpha0 + Np/2, beta0 + SSBt/2)
# sample the parameter fropm posterior
theta[iter, t] = sliceSample(n = 1, target, x.interval=c(trt, upperBound), root.accuracy = upperBound/1000)







# parameters to monitor
Niter <- 100
sigma <- matrix(nrow=Niter, ncol=nItem)
SigmaFull <- array(dim=c(Niter, Nt+1, nItem, nItem))
theta <- matrix(nrow=Niter, ncol=1)

# initial values
sigma[1,] <- runif(nItem, 0.5, 1.5) # item unique variances
SigmaFull[1,1,,] <- diag(sigma[1,], ncol=nItem, nrow=nItem) # full covariance matrix at different layers
theta[1,] <- 0 # covariance parameter


# priors
alpha0 <- 0.1
beta0 <- 0.1
# iterators
iter <- i <- t <- 1
# loop
for(iter in 1:Niter){
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
  SSWi <- SSW[i,i]
  # ====================
  #   Sample sigma[i]
  Q <- theta[iter,t]*u[i,t]
  target <- function(x){ dIG(x, alpha0 + Np/2, beta0 + SSWi/2, Q, 0) }
  upperBound <- 1/qgamma(.001, alpha0 + Np/2, beta0 + SSWi/2)
  sigma[iter, i] = sliceSample(n = 1, target, x.interval=c(0, upperBound), root.accuracy = upperBound/1000)
  # Computer SSB eq (8)
  Yut <- Y[,which(u[,t] != 0)]
  Y.t <- matrix(rep(mean(Yut), nrow(Yut)), ncol=1)
  Y.it <- matrix(rowMeans(Yut), ncol=1)
  SSBt <- c(t(Y.it - Y.t)%*%(Y.it - Y.t))
  # =====================
  # shift parameter
  psi <- SSBt - theta[iter, t]
  # =====================
  # truncation parameter
  SigmaFull[iter,1, , ] <- diag(sigma[iter, ], ncol=nItem, nrow=nItem)
  SigmaFull[iter, t+1, , ] <- diag(sigma[iter, ], ncol=nItem, nrow=nItem) + theta[iter, t] * (u[,t, drop=F]%*%t(u[,t, drop=F]))
  # get iterative inverses
  SFt <- SigmaFull[iter,t, , ]
  SFtp1 <- SigmaFull[iter,t+1, , ]
  # inverse of full matrix
  SF0inv <- solve(SigmaFull[iter, 1, , ])
  # compute denominator of eq (5)
  invtheta <- ifelse(theta[t] == 0, 0, 1/theta[t])
  devisor <- c(invtheta + t(u[,t, drop=F])%*%SF0inv%*%u[,t, drop=F])
  # compute Inverse at layer t
  SFtinv <- SF0inv - (SF0inv%*%u[,t, drop=F]%*%t(u[,t, drop=F])%*%SF0inv)/( devisor )

  # compute lower bound of theta[t]
  trt <- -1/(t(u[,t, drop=F])%*%SFtinv%*%u[,t, drop=F])
  if(trt > 0) trt <- 0

  # =====================
  #   Sample theta
  # set up posterior
  target <- function(x){ dIG(x, alpha0 + Np/2, beta0 + SSBt/2, psi, trt) }
  # compute an upper bound for slice sampling
  upperBound <- 1/qgamma(.001, alpha0 + Np/2, beta0 + SSBt/2)
  # sample the parameter fropm posterior
  theta[iter, t] = sliceSample(n = 1, target, x.interval=c(trt, upperBound), root.accuracy = upperBound/1000)


}

