source("code/load_packages.R")
options(max.print = 10000, scipen = 10, digits=2)

# useful functions
invlogit <- function(x) {exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
# Generating Data
N <- 200 # number of respondents
J <- 5 # number of items
C <- 3 # number of response categories
# ========================= #
# latent person parameters
etaCor <- 0.23 # correlation between ability and speediness
etasd <- c(1,sqrt(0.1))
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

# set up estimation
mydata <- list(
  y = as.matrix(Y),
  N = nrow(Y),
  J=5, C=3
)
jags.model <- function(){
  ### Model
  for(i in 1:N){
    for(j in 1:J){
      # data model
      y[i,j] ~ dcat(omega[i,j, ])

      # LRV
      ystar[i,j] ~ dnorm(lambda[j]*eta[i], 1)

      # Pr(nu = 3)
      pi[i,j,3] = phi(ystar[i,j] - tau[j,2])
      # Pr(nu = 2)
      pi[i,j,2] = phi(ystar[i,j] - tau[j,1]) - phi(ystar[i,j] - tau[j,2])
      # Pr(nu = 1)
      pi[i,j,1] = 1 - phi(ystar[i,j] - tau[j,1])

      # observed category prob (Pr(y=c))
      for(c in 1:C){
        omega[i,j, c] = gamma[1,c,j]*pi[i,j,1] +
          gamma[2,c,j]*pi[i,j,2] +
          gamma[3,c,j]*pi[i,j,3]
      }
    }
  }
  ### Misclassification parameters
  for(j in 1:J){
    gamma[1,1:3,j] ~ ddirch(alpha[1,1:3])
    gamma[2,1:3,j] ~ ddirch(alpha[2,1:3])
    gamma[3,1:3,j] ~ ddirch(alpha[3,1:3])
  }
  ### Priors
  # person parameters
  for(i in 1:N){
    eta[i] ~ dnorm(0, 1)
  }
  # item parameters
  for(j in 1:J){
    # Thresholds
    tau[j, 1] ~ dnorm(0.0,0.1)
    tau[j, 2] ~ dnorm(0, 0.1);T(tau[j, 1],)
    # loadings
    lambda[j] ~ dbeta(2,2)
    ### LRV
    theta[j] = 1 + pow(lambda[j],2)
  }
}

# vector of all parameters to save
#param_save <- c("psi", "lambda", "tau","V", "p", "ystar")
param_save <- c("lambda",  "tau","theta", "gamma", "pi", "omega", "yppc")

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = param_save,
  n.iter=2000,
  #n.thin = 50,
  #n.burnin = 0000,
  n.chains = 4
)

#print(fit)

fit$BUGSoutput$summary[paste0("lambda[",1:5,"]") ,]
fit$BUGSoutput$summary[paste0("theta[",1:5,"]") ,]
fit$BUGSoutput$summary[c(paste0("tau[1,",1:2,"]"), paste0("tau[2,",1:2,"]"), paste0("tau[3,",1:2,"]"), paste0("tau[4,",1:2,"]"), paste0("tau[5,",1:2,"]")),]

dat[1,]
fit$BUGSoutput$summary[paste0("pi[1,1,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("pi[1,2,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("pi[1,3,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("pi[1,4,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[1,1,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[1,2,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[1,3,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[1,4,",1:3,"]") ,]

dat[42,]
fit$BUGSoutput$summary[paste0("pi[42,1,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("pi[42,2,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("pi[42,3,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("pi[42,4,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[42,1,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[42,2,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[42,3,",1:3,"]") ,]
fit$BUGSoutput$summary[paste0("omega[42,4,",1:3,"]") ,]


# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(
  plot.data,
  pars = c(paste0("tau[1,",1:2,"]"),
           paste0("tau[2,",1:2,"]"),
           paste0("tau[3,",1:2,"]"),
           paste0("tau[4,",1:2,"]"),
           paste0("tau[5,",1:2,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("lambda[",1:5,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("theta[",1:5,"]")),
  prob = 0.8) +
  plot_title


mcmc_areas(
  plot.data,
  pars = c(paste0("gamma[1,1," ,1:5,"]"),
           paste0("gamma[1,2," ,1:5,"]"),
           paste0("gamma[1,3," ,1:5,"]"),
           paste0("gamma[2,1," ,1:5,"]"),
           paste0("gamma[2,2," ,1:5,"]"),
           paste0("gamma[2,3," ,1:5,"]"),
           paste0("gamma[3,1," ,1:5,"]"),
           paste0("gamma[3,2," ,1:5,"]"),
           paste0("gamma[3,3," ,1:5,"]")
  ),
  prob = 0.8) +
  plot_title




# set up estimation - using the code from Molenaar et al. (2021)
mydata <- list(
  x = Y,
  lrt = logt,
  N = nrow(Y),
  nit=J, ncat=C
)
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # GRM model
      for(k in 2:(ncat)){
        # P(greater than or equal to category k > 1)
        P.gte[p,i,k] <- exp(a[i]*theta[p]+b[i,k])/(1+exp(a[i]*theta[p]+b[i,k]))
      }
      # P(greater than or equal to category 1)
      P.gte[p,i,1] <- 1
      # equal to prob.
      for(k in 1:(ncat-1)){
        # P(greater equal to category k < K)
        Prob[p,i,k] <- P.gte[p,i,k]-P.gte[p,i,k+1]
      }
      # P(greater equal to category K)
      Prob[p,i,ncat] <- P.gte[p,i,ncat]
      # observable dist.
      x[p,i] ~ dcat(Prob[p,i,1:ncat])

      # log-RT model
      dev[p,i]<-a[i]*(theta[p] - (b[i,2]+b[i,3]+b[i,4])/3)
      mu.rt[p,i]<- icept[i] - speed[p] - rho2*dev[p,i] - rho1 * abs(dev[p,i]-delta)
      lrt[p,i]~dnorm(mu.rt[p,i],prec[i])
    }}
  for(p in 1:N){
    theta[p]~dnorm(0,1)
    speed[p]~dnorm(sigma.ts*theta[p],prec.s)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # item disc.
    a[i]~dnorm(0,.1); I(0,)
    # location parameters
    b[i,2] ~ dnorm(2, .5)
    b[i,3] ~ dnorm(1, .5);I( , b[i,2])
    b[i,4] ~ dnorm(1, .5);I( , b[i,3])
  }
  negrho <- -1*rho1
  rho1~dnorm(0,.1);I(0,)
  rho2~dnorm(0,.1);I(negrho,rho1)
  delta~dnorm(0,.1)

  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))

} # end jags model



jags.params <- c("delta", "rho1", "rho2", "b", "a", "icept","prec.s", "sigma.t", "sigma.ts", "cor.ts")


# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = jags.params,
  n.iter=2000,
  #n.thin = 50,
  #n.burnin = 0000,
  n.chains = 4
)
fit1 <- fit
fit <- fit1
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(
  plot.data,
  pars = c(paste0("b[1,",2:C,"]"),
           paste0("b[2,",2:C,"]"),
           paste0("b[3,",2:C,"]"),
           paste0("b[4,",2:C,"]"),
           paste0("b[5,",2:C,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("a[",1:J,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c("sigma.t", "sigma.ts", "cor.ts"),
  prob = 0.8) +
  plot_title



# misclassification model
jags.model <- function(){
  ### Model
  for(p in 1:N){
    for(i in 1:nit){
      # GRM model
      for(k in 2:(ncat)){
        # P(greater than or equal to category k > 1)
        P.gte[p,i,k] <- exp(a[i]*theta[p]+b[i,k])/(1+exp(a[i]*theta[p]+b[i,k]))
      }
      # P(greater than or equal to category 1)
      P.gte[p,i,1] <- 1
      # equal to prob.
      for(k in 1:(ncat-1)){
        # P(greater equal to category k < K)
        Prob[p,i,k] <- P.gte[p,i,k]-P.gte[p,i,k+1]
      }
      # P(greater equal to category K)
      Prob[p,i,ncat] <- P.gte[p,i,ncat]
      # observable dist.
      x[p,i] ~ dcat(omega[p,i,1:ncat])


      # log-RT model
      dev[p,i]<-a[i]*(theta[p] - (b[i,2]+b[i,3]+b[i,4])/3)
      mu.rt[p,i]<- icept[i] - speed[p] - rho2*dev[p,i] - rho1 * abs(dev[p,i]-delta)
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
        omega[p,i, c] = gamma[p,i,1,c]*Prob[p,i,1] + gamma[p,i,2,c]*Prob[p,i,2] + gamma[p,i,3,c]*Prob[p,i,3]
      }

    }
  }

  # person parameters
  for(p in 1:N){
    theta[p]~dnorm(0,1)
    speed[p]~dnorm(sigma.ts*theta[p],prec.s)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    prec[i]~dgamma(.1,.1)
    # item disc.
    a[i]~dnorm(0,.1); I(0,)
    # location parameters
    b[i,2] ~ dnorm(2, .5)
    b[i,3] ~ dnorm(1, .5);I( , b[i,2])
    b[i,4] ~ dnorm(1, .5);I( , b[i,3])
  }
  negrho <- -1*rho1
  rho1~dnorm(0,.1);I(0,)
  rho2~dnorm(0,.1);I(negrho,rho1)
  delta~dnorm(0,.1)

  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))

} # end jags model



# vector of all parameters to save
#param_save <- c("psi", "lambda", "tau","V", "p", "ystar")
# param_save <- c("lambda",  "tau","theta", "gamma", "pi", "omega",
#                 "delta", "rho1", "rho2",  "icept","prec.s",
#                 "sigma.t", "sigma.ts", "cor.ts",
#                 "eta", "speed")
# param_save <- c("lambda",  "tau","theta", "gamma",
#                 "delta", "rho1", "rho2",  "icept","prec.s",
#                 "sigma.t", "sigma.ts", "cor.ts")
jags.params <- c("delta", "rho1", "rho2", "b", "a", "icept","prec.s", "sigma.t", "sigma.ts", "cor.ts", "gamma")
# data
mydata <- list(
  x = Y,
  lrt = logt,
  N = nrow(Y),
  nit=J, ncat=C
)

# fit model
fit <- jags(
  model.file=jags.model,
  data=mydata,
  #inits = start,
  parameters.to.save = jags.params,
  n.iter=2000,
  #n.thin = 50,
  #n.burnin = 0000,
  n.chains = 4
)
fit2 <- fit
print(fit)
fit <- fit2
fit$BUGSoutput$summary[! rownames(fit$BUGSoutput$summary) %like% "gamma",]

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)


plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(
  plot.data,
  pars = c(paste0("b[1,",2:C,"]"),
           paste0("b[2,",2:C,"]"),
           paste0("b[3,",2:C,"]"),
           paste0("b[4,",2:C,"]"),
           paste0("b[5,",2:C,"]")),
  prob = 0.8) +
  plot_title

mcmc_areas(
  plot.data,
  pars = c(paste0("a[",1:J,"]")),
  prob = 0.8) +
  plot_title+
  lims(x=c(0, 3))

mcmc_areas(
  plot.data,
  pars = c("sigma.t", "sigma.ts", "cor.ts"),
  prob = 0.8) +
  plot_title



# misc. test
jags.model <- function(){
  ### Model
  for(i in 1:N){
    for(j in 1:J){
      # data model
      y[i,j] ~ dcat(omega[i,j, ])

      # LRV
      ystar[i,j] ~ dnorm(lambda[j]*eta[i], invtheta[j])

      # Pr(nu = 3)
      pi[i,j,3] = phi(ystar[i,j] - tau[j,2])
      # Pr(nu = 2)
      pi[i,j,2] = phi(ystar[i,j] - tau[j,1]) - phi(ystar[i,j] - tau[j,2])
      # Pr(nu = 1)
      pi[i,j,1] = 1 - phi(ystar[i,j] - tau[j,1])

      # observed category prob (Pr(y=c))
      for(c in 1:C){
        omega[i,j, c] = gamma[i,j,1,c]*pi[i,j,1] + gamma[i,j,2,c]*pi[i,j,2] + gamma[i,j,3,c]*pi[i,j,3]
      }
      # log-RT model
      dev[i,j]<- eta[i] - (tau[j,1]+tau[j,2])/2
      mu.rt[i,j]<- icept[j] - speed[i] - rho2*dev[i,j] - rho1 * abs(dev[i,j]-delta)
      lrt[i,j]~dnorm(mu.rt[i,j],prec[j])

      # misclassification parameters (as fixed functions of latent RT)
      for(c in 1:C){
        for(ct in 1:C){
          gamma[i,j,c,ct] <- ifelse(c == ct,
                                    ilogit(icept[j] - speed[i]),
                                    (1/(C-1))*(1-ilogit(icept[j] - speed[i]))
          )
        }
      }

    }
  }
  ### LRV
  for(j in 1:J){
    theta[j] = 1 - pow(lambda[j],2)
    invtheta[j] = pow(theta[j], -1)
  }
  ### Priors
  # person parameters
  for(i in 1:N){
    eta[i] ~ dnorm(0, 1)
    speed[i]~dnorm(sigma.ts*eta[i], prec.s)
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s~dgamma(.1,.1)
  for(j in 1:J){
    # lrt parameters
    icept[j]~dnorm(0,.1)
    prec[j]~dgamma(.1,.1)
    # Thresholds
    tau[j, 1] ~ dnorm(0.0,0.1)
    tau[j, 2] ~ dnorm(0, 0.1);T(tau[j, 1],)
    # loadings
    lambda[j] ~ dbeta(2,2)
  }
  # response time parameters
  negrho <- -1*rho1
  rho1~dnorm(0,.1);I(0,)
  rho2~dnorm(0,.1);I(negrho,rho1)
  delta~dnorm(0,.1)

  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))
}
