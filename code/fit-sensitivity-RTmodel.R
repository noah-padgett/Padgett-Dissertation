source("code/load_packages.R")
options(max.print = 10000, scipen = 10, digits=2)

# useful functions
invlogit <- function(x) {exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
# Generating Data
N <- 200 # number of respondents
J <- 5 # number of items
C <- 2 # number of response categories
# ========================= #
# latent person parameters
etaCor <- 0.5 # correlation between ability and speediness
etasd <- c(1,sqrt(0.1))
eta <- mvtnorm::rmvnorm(N, mean = c(0, 0), sigma = matrix(c(etasd[1], etasd[2]*etaCor, etasd[2]*etaCor, etasd[2]**2), ncol = 2))
eta0 <- matrix(eta[,1],nrow=1) # ability
eta1 <- matrix(eta[,2],nrow=1) # log speediness
# ========================= #
# item parameters
# item factor loadings
lambda <- matrix(rep(0.9, J), ncol=1)
# item latent residual variances
theta <- c(1 - lambda**2)
# item thresholds
tau <- matrix(ncol=C-1, nrow=J)
for(c in 1:(C-1)){
  if(c == 1){
    tau[,1] <- runif(J, -1, -0.33)
  }
  if(c > 1){
    tau[,c] <- tau[,c-1] + runif(J, 0.25, 1)
  }
}
# latent item response
ystar <- lambda%*%eta0
ystar <- apply(ystar, 2, FUN = function(x){mvtnorm::rmvnorm(1, x, diag(theta, ncol=J, nrow=J))})
# response time parameters (copied from Molenaar et al. 2021)
nu <- matrix(rep(2, J), ncol=1)
sigma.ei <- matrix(rep(0.25, J), ncol=1)
rho1 <- 0.1
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

for(i in 1:N){
  for(j in 1:J){
  for(b in 1:C){
  for(c in 1:C){
    gamma[i,j,b,c] <- misclass.time.trans(nu[j, 1] - eta1[1, i],
                                          b, c, C)
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
      #omega[i, j, c] = gamma[i,j,1,c]*pi[i,j,1] + gamma[i,j,2,c]*pi[i,j,2] + gamma[i,j,3,c]*pi[i,j,3]
    }
    Y[i,j] <- sample(x=1:C, size=1, replace=F, prob=omega[i,j,])
  }
}

# set up estimation
mydata <- list(
  x = Y,
  lrt = logt,
  N = nrow(Y),
  nit=5, ncat=2
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # Rasch model
      logit(Prob[p,i,1]) <- theta[p]-b[i]
      Prob[p,i,2] <- 1 - Prob[p,i,1]
      # observable dist.
      x[p,i] ~ dbern(omega[p,i,1])

      # log-RT model
      dev[p,i]<-theta[p] - b[i]
      mu.rt[p,i]<- icept[i] - speed[p] - abs(dev[p,i])
      lrt[p,i]~dnorm(mu.rt[p,i],prec[i])

      # misclassificaiton rates
      for(c in 1:ncat){
        for(ct in 1:ncat){
          gamma[p,i,c,ct] <- ifelse(c == ct,
                                    ilogit(icept[i] - speed[p]),
                                    (1/(ncat-1))*(1-ilogit(icept[i] - speed[p]))
          )
        }
      }

      # misclassification weighting
      for(c in 1:ncat){
        omega[p,i, c] = gamma[p,i,1,c]*Prob[p,i,1] + gamma[p,i,2,c]*Prob[p,i,2]
      }

    }
  }
  # person parameters
  for(p in 1:N){
    theta[p]~dlogis(0,1)
    speed[p]~dnorm(sigma.ts*theta[p],prec.s)
  }

  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # location parameters
    b[i] ~ dnorm(0, .5)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)

  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))
  # Rasch item fit statistics
  # deviations
  for(i in 1:nit){
    for(p in 1:N){
      w[p,i] <- omega[p,i,1]*(1-omega[p,i,1])
      z[p,i] <- (x[p,i]-omega[p,i,1])/pow(w[p,i], 0.5)
      wz2[p,i] <- w[p,i]*pow(z[p,i],2)
    }
    outfitMSQ[i] <- sum(pow(z[,i], 2))/N
    infitMSQ[i] <- sum(wz2[,i])/sum(w[,i])
  }

} # end jags model



# vector of all parameters to save
jags.params <- c("b", "icept","prec","prec.s", "sigma.t", "sigma.ts", "cor.ts", "outfitMSQ", "infitMSQ")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=10000,
  n.chains = 4
)
fit1 <- fit
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(
  plot.data,
  pars = c(paste0("b[",1:5,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("outfitMSQ[",1:5,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("infitMSQ[",1:5,"]")),
  prob = 0.8) +
  plot_title


# reduced model without RT
mydata <- list(
  x = Y,
  N = nrow(Y),
  nit=5
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # Rasch model
      logit(Prob[p,i]) <- theta[p]-b[i]
      # observable dist.
      x[p,i] ~ dbern(Prob[p,i])

    }
  }
  # person parameters
  for(p in 1:N){
    theta[p]~dlogis(0,1)
  }

  for(i in 1:nit){
    # location parameters
    b[i] ~ dnorm(0, .5)
  }
  # Rasch item fit statistics
  # deviations
  for(i in 1:nit){
    for(p in 1:N){
      w[p,i] <- Prob[p,i]*(1-Prob[p,i])
      z[p,i] <- (x[p,i]-Prob[p,i])/pow(w[p,i], 0.5)
      wz2[p,i] <- w[p,i]*pow(z[p,i],2)
    }
    outfitMSQ[i] <- sum(pow(z[,i], 2))/N
    infitMSQ[i] <- sum(wz2[,i])/sum(w[,i])
  }

} # end jags model



# vector of all parameters to save
jags.params <- c("b", "outfitMSQ", "infitMSQ")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  parameters.to.save = jags.params,
  n.iter=10000,
  n.chains = 4
)
fit2 <- fit
print(fit)
jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(
  plot.data,
  pars = c(paste0("b[",1:5,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("outfitMSQ[",1:5,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("infitMSQ[",1:5,"]")),
  prob = 0.8) +
  plot_title



# eRm
library(eRm)
fit.R <- eRm::RM(Y)
summary(fit.R)
pp <- person.parameter(fit.R)
itemfit(pp)
