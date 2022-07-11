# covariance structure modeling
# test
source("code/load_packages.R")
options(max.print = 10000, scipen = 10, digits=2)

mydata <- readr::read_csv("data/person.csv",col_names = F)
colnames(mydata) <- c(paste0("q",0:49),"gender","age","E","C","N","O","A")
# item grps
itemExt = c("q0", "q5", "q10", "q15", "q20", "q25", "q30", "q35", "q40", "q45")
itemNeu = c("q1", "q6", "q11", "q16", "q21", "q26", "q31", "q36", "q41", "q46")
itemAgr = c("q2", "q7", "q12", "q17", "q22", "q27", "q32", "q37", "q42", "q47")
itemCon = c("q3", "q8", "q13", "q18", "q23", "q28", "q33", "q38", "q43", "q48")
itemOpn = c("q4", "q9", "q14", "q19", "q24", "q29", "q34", "q39", "q44", "q49")

# 0's are NA so remove
mydata[mydata == 0] <- NA
mydata <- na.omit(mydata)
# random subset
set.seed(1)
ids <- sample(1:nrow(mydata), 100)
# Prep Data - extraversion
dat <- mydata[ids,itemExt] # 5 categories (1, 2, 3, 4, 5)

a <- dat %>%
  summarise(across(everything(),.fns = mean))
a[2,] <- dat %>%
  summarise(across(everything(),.fns = sd))
a[3,] <- dat %>%
  summarise(across(everything(),.fns = min))
a[4,] <- dat %>%
  summarise(across(everything(),.fns = max))
a <- a %>%  as.matrix()
rownames(a) <- c("mean", "sd", "min", "max")
t(a)

# recode
reVar <- c("q5", "q15","q25", "q35", "q45")
#, "q6", "q11", "q16", "q21", "q26", "q31", "q36", "q41", "q46","q7", "q12", "q17", "q22", "q27", "q32", "q37", "q42", "q47", "q8", "q13", "q18", "q23", "q28", "q33", "q38", "q43", "q48", "q9", "q14", "q19", "q24", "q29", "q34", "q39", "q44", "q49")
dat2 <- dat
dat2[,reVar] <- apply(
  dat[,reVar],
  2,
  FUN = function(x){
    as.numeric(fct_rev(as.factor(x)))
  }
)

apply(dat, 2, table)
apply(dat2, 2, table)

nj <- 10
m <- 100
dat$id <- 1:m
datlong <- dat %>% pivot_longer(cols=contains("q"), names_to="item", values_to="y")
mydata <- list(
  nj = nj,
  m = m,
  id = matrix(1:nrow(datlong), ncol=nj, nrow=m, byrow=T),
  Imat = diag(nrow=nj, ncol=nj),
  Ivec = matrix(1, nrow=nj, ncol=1),
  Jmat = matrix(1, ncol=nj, nrow=nj),
  Y = as.matrix(datlong[,3]) #, X = dat[,2,drop=F]
)

# model
# model
bcsm.model <- function(){

  for(j in 1:m){
    omega[j, 1:nj , 1:nj] <- Imat*sigma2 + Jmat*tau[j]
    beta.vec[j,1:nj] <- Ivec*beta[1]
    Y[id[j,],] ~ dmnorm.vcov(beta.vec[j,], omega[j,,])
  }
  # parameters
  #inv.sigma ~ dunif(0,3) #dgamma(5, 10) # error precision
  sigma2 ~ dunif(0,3) #<- 1/inv.sigma      # reparameterize to variance

  lb <- -sigma2/nj + 0.001
  tau.mu ~ dnorm(0, 0.5);T(lb,) # truncated normal for  covariance
  tau.sig ~ dt(0, 1, nj);T(0,) #dgamma(1, 20) #dunif(0,0.2) #
  tau.prec = pow(tau.sig, -2)
  for(j in 1:m){
    tau[j] ~ dnorm(tau.mu, tau.prec);T(lb,)
  }
  beta[1] ~ dnorm(0, 0.1)

}
# initial values
start_values <- list(
  list("sigma2"=1,
       "tau.mu"=0,
       "tau.sig"=0.05,
       "tau"=rep(0,m),
       "beta"=rep(0, nj))
)

# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("sigma2", "tau.mu", "tau.sig", "tau", "beta")

# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  #inits=start_values,
  parameters.to.save = param_save,
  n.iter=10000,n.chains = 4)

print(fit0)
fit <- fit0

fit <- R2jags::autojags(fit0,n.iter = 50000, n.thin = 20, Rhat = 1.01)
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,regex_pars = "tau.mu", prob=0.8)
mcmc_areas(plot.data,regex_pars = "beta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma2", prob=0.8)



bcsm.model <- function(){

  for(j in 1:m){
    omega[j, 1:nj , 1:nj] <- Imat*sigma2 + Jmat*tau
    beta.vec[j,1:nj] <- Ivec*beta[1]
    Y[id[j,],] ~ dmnorm.vcov(beta.vec[j,], omega[j,,])
  }
  # parameters
  #inv.sigma ~ dunif(0,3) #dgamma(5, 10) # error precision
  sigma2 ~ dunif(0,3) #<- 1/inv.sigma      # reparameterize to variance
  lb <- -sigma2/nj + 0.001
  tau ~ dnorm(0, 0.5);T(lb,) # truncated normal for  covariance
  #for(i in 1:nj){
    beta[1] ~ dnorm(0, 0.1)
  #}


}
# initial values
start_values <- list(
  list("sigma2"=1,
       "tau"=0,
       "beta"=0 )#rep(0, nj))
)

# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("sigma2", "tau", "beta")

# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  #inits=start_values,
  parameters.to.save = param_save,
  n.iter=1000,n.chains = 2)

print(fit0)
fit <- fit0

fit <- R2jags::autojags(fit0,n.iter = 50000, n.thin = 20, Rhat = 1.01)
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,regex_pars = "tau", prob=0.8)
mcmc_areas(plot.data,regex_pars = "beta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma2", prob=0.8)



# Matrix version
itemExt = c("q0", "q5", "q10", "q15", "q20", "q25", "q30", "q35", "q40", "q45")

u <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), nrow=10, ncol=1)
diagSigma <- matrix(0, ncol=nj, nrow=nj)
diag(diagSigma) <- NA
mydata <- list(
  nj = nj,
  m = m,
  id = matrix(1:nrow(datlong), ncol=nj, nrow=m, byrow=T),
  u=u,
  Y = as.matrix(datlong[,3]),
  diagSigma = diagSigma
)

bcsm.model <- function(){

  for(j in 1:m){
    Y[id[j,],] ~ dmnorm.vcov(beta[], omega[,])
  }
  # error variance
  for(i in 1:nj){
    sigma[i] ~ dunif(0,3)
    diagSigma[i,i] <- sigma[i]
  }
  SigmaInv[1:nj, 1:nj] <- inverse(diagSigma[1:nj, 1:nj] )
  # obtain lower bound - truncation
  tr[1] <- - 1 /(t(u) %*% SigmaInv %*% u)
  # sample thet
  theta[1] ~ dnorm(0, 0.5);T(tr[1],)
  # compute u vector products
  ut <- inprod(u,t(u))
  # update full covariance matrix
  omega[1:nj , 1:nj] <- diagSigma + ut*theta[1]
  # item intercepts
  for(i in 1:nj){
    beta[i] ~ dnorm(0, 0.1)
  }




  # for(t in 2:(Nt+1)){
  #   # check value of theta for approx 0 value
  #   LAM <- ifelse(abs(theta[t-1]) < 1e-5, 1e5, 1/theta[t-1])
  #   d <- c(LAM + inprod((inprod(t(u[,t-1]),SigmaInv[t-1, 1:nj, 1:nj]),u[,t-1]))
  #   # update layer (t) inverse
  #   SigmaInv[t, 1:nj, 1:nj] <- SigmaInv[t-1, 1:nj, 1:nj] - (inprod(inprod(inprod(SigmaInv[t-1, 1:nj, 1:nj],u[,t-1]),t(u[,t-1])),SigmaInv[t-1, 1:nj, 1:nj]))/d
  #
  #   # update truncation parameter
  #   tr[t-1] <- - 1/(inprod(inprod(t(u[,t-1]),SigmaInv[t-1, 1:nj, 1:nj]),u[,t-1]))
  # }

  # # sample theta[t]
  # for(t in 1:Nt){
  #   theta[t] ~ dnorm(0, 0.5);T(tr[t],)
  # }

}

# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("sigma","beta","theta")

# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  #inits=start_values,
  parameters.to.save = param_save,
  n.iter=1000,n.chains = 1)

print(fit0)
fit <- fit0

fit <- R2jags::autojags(fit0,n.iter = 50000, n.thin = 20, Rhat = 1.01)
print(fit)

jags.mcmc <- as.mcmc(fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,regex_pars = "tau", prob=0.8)
mcmc_areas(plot.data,regex_pars = "beta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma2", prob=0.8)


# model with factor loadings in u
u <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              0, 1, 0, 1, 0, 1, 0, 1, 0, 1), ncol=10, nrow=2, byrow=T)

# expanded u matrix
itemExt = c("q0", "q5", "q10", "q15", "q20", "q25", "q30", "q35", "q40", "q45")
reVar <- c("q5", "q15","q25", "q35", "q45")

u <- matrix(c(rep(0,nj),
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              0, 1, 0, 1, 0, 1, 0, 1, 0, 1), nrow=10, ncol=3)
diagSigma <- matrix(0, ncol=nj, nrow=nj)
diag(diagSigma) <- NA
mydata <- list(
  nj = nj,
  m = m,
  id = matrix(1:nrow(datlong), ncol=nj, nrow=m, byrow=T),
  u=u,
  Nt = 2,
  Y = as.matrix(datlong[,3]),
  diagSigma = diagSigma
)

bcsm.model <- function(){

  for(j in 1:m){
    Y[id[j,],] ~ dmnorm.vcov(beta[], omega[,])
  }
  # error variance
  for(i in 1:nj){
    sigma[i] ~ dunif(0,3)
    diagSigma[i,i] <- sigma[i]
  }
  SigmaInv[1, 1:nj, 1:nj] <- inverse(diagSigma[1:nj, 1:nj] )

  # =========================
  # layers
  # check value of theta for approx 0 value
  for(t in 2:(Nt+1)){
    LAM[t] <- ifelse(abs(theta[t-1]) < 1e-5, 1e5, 1/theta[t-1])
    d[t] <- LAM[t] + t(u[,t])%*%SigmaInv[t-1,1:nj, 1:nj]%*%u[,t]
    # update layer (t) inverse
    A[1:nj, t] <- SigmaInv[t-1,1:nj, 1:nj]%*%u[1:nj,t]
    for(i in 1:nj){
      for(j in 1:nj){
        B[t,i,j] <- A[i,t]*u[j,t]
      }
    }

    C[t,1:nj,1:nj] <- B[t, , ]%*%SigmaInv[t-1, , ]
    SigmaInv[t, 1:nj, 1:nj] <- SigmaInv[t-1, 1:nj, 1:nj] - (C[t,,])/d[t]
    # update truncation parameter
    tr[t-1] <- - 1/(t(u[,t])%*%SigmaInv[t-1, 1:nj, 1:nj]%*%u[,t])
  }
  # # sample theta[t]
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
  # update full covariance matrix
  omega[1:nj , 1:nj] <- diagSigma + utsum

  # item intercepts
  for(i in 1:nj){
    beta[i] ~ dnorm(0, 0.1)
  }

}


# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("sigma","beta","theta", "tr")

# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  #inits=start_values,
  parameters.to.save = param_save,
  n.iter=1000,n.chains = 4)

print(fit0)
fit <- fit0

fit <- R2jags::autojags(fit0,n.iter = 50000, n.thin = 20, Rhat = 1.01)
print(fit)

jags.mcmc <- as.mcmc(fit)
#jags.mcmc <- as.mcmc(fit$BUGSoutput$sims.matrix)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,pars = "theta[1]", prob=0.8)
mcmc_areas(plot.data,pars = "theta[2]", prob=0.8)
mcmc_areas(plot.data,regex_pars = "beta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma", prob=0.8)

mcmc_acf(plot.data, regex_pars = "theta")
mcmc_acf(plot.data, regex_pars = "beta")
mcmc_acf(plot.data, regex_pars = "sigma")

mcmc_pairs(plot.data, regex_pars = "theta")

mcmc_trace(plot.data, regex_pars = "theta")
mcmc_trace(plot.data, regex_pars = "beta")
mcmc_trace(plot.data, regex_pars = "sigma")

gelman.plot(jags.mcmc)




# Non-unit factor factor loadings for main factor
itemExt = c("q0", "q5", "q10", "q15", "q20", "q25", "q30", "q35", "q40", "q45")
reVar <- c("q5", "q15","q25", "q35", "q45")

u <- matrix(c(rep(0,nj),
              rep(1,10),
              0, 1, 0, 1, 0, 1, 0, 1, 0, 1), nrow=10, ncol=3)
diagSigma <- matrix(0, ncol=nj, nrow=nj)
diag(diagSigma) <- NA
mydata <- list(
  nj = nj,
  m = m,
  id = matrix(1:nrow(datlong), ncol=nj, nrow=m, byrow=T),
  #u=u,
  Nt = 1,
  Y = as.matrix(datlong[,3]),
  diagSigma = diagSigma
)

bcsm.model <- function(){

  for(j in 1:m){
    # data model
    Y[id[j,],] ~ dmnorm.vcov(beta[], omega[,])
  }
  # factor loadings
  for(i in 1:nj){
    lambda[i] ~ dnorm(0, 1)
  }
  # item intercepts
  for(i in 1:nj){
    beta[i] ~ dnorm(0, 0.1)
  }

  # error variance
  for(i in 1:nj){
    sigma[i] ~ dunif(0,3)
    diagSigma[i,i] <- sigma[i]
  }
  SigmaInv[1, 1:nj, 1:nj] <- inverse(diagSigma[1:nj, 1:nj] )

  # compute truncate para
  u[1:nj,1] <- lambda[]
  tr[1] <- - 1/(t(u[1:nj,1])%*%SigmaInv[1, 1:nj, 1:nj]%*%u[1:nj,1])

  # =========================
  # layers
  # check value of theta for approx 0 value
  # for(t in 2:(Nt+1)){
  #   LAM[t] <- ifelse(abs(theta[t-1]) < 1e-5, 1e5, 1/theta[t-1])
  #   d[t] <- LAM[t] + t(u[,t])%*%SigmaInv[t-1,1:nj, 1:nj]%*%u[,t]
  #   # update layer (t) inverse
  #   A[1:nj, t] <- SigmaInv[t-1,1:nj, 1:nj]%*%u[1:nj,t]
  #   for(i in 1:nj){
  #     for(j in 1:nj){
  #       B[t,i,j] <- A[i,t]*u[j,t]
  #     }
  #   }
  #
  #   C[t,1:nj,1:nj] <- B[t, , ]%*%SigmaInv[t-1, , ]
  #   SigmaInv[t, 1:nj, 1:nj] <- SigmaInv[t-1, 1:nj, 1:nj] - (C[t,,])/d[t]
  #   # update truncation parameter
  #   tr[t-1] <- - 1/(t(u[,t])%*%SigmaInv[t-1, 1:nj, 1:nj]%*%u[,t])
  # }
  # # sample theta[t]
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
  # update full covariance matrix
  omega[1:nj , 1:nj] <- diagSigma + utsum

}


# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("lambda", "sigma","beta","theta", "tr")

# helps to start all factor loadings as 1s
# makes the first iteration estimable
start_values <- list(
  list("lambda"=rep(1, nj)),
  list("lambda"=rep(1, nj))
)

# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  inits=start_values,
  parameters.to.save = param_save,
  n.iter=1000,n.chains = 2)

print(fit0)
fit <- fit0

fit <- R2jags::autojags(fit0,n.iter = 50000, n.thin = 20, Rhat = 1.01)
print(fit)

jags.mcmc <- as.mcmc(fit)
#jags.mcmc <- as.mcmc(fit$BUGSoutput$sims.matrix)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,pars = "theta[1]", prob=0.8)
mcmc_areas(plot.data,pars = "theta[2]", prob=0.8)
mcmc_areas(plot.data,regex_pars = "beta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma", prob=0.8)

mcmc_acf(plot.data, regex_pars = "theta")
mcmc_acf(plot.data, regex_pars = "beta")
mcmc_acf(plot.data, regex_pars = "sigma")

mcmc_pairs(plot.data, regex_pars = "theta")

mcmc_trace(plot.data, regex_pars = "theta")
mcmc_trace(plot.data, regex_pars = "beta")
mcmc_trace(plot.data, regex_pars = "sigma")

gelman.plot(jags.mcmc)


# Non-unit factor factor loadings for main factor
u <- matrix(c(rep(0,nj),
              rep(NA,10),
              0, 1, 0, 1, 0, 1, 0, 1, 0, 1), nrow=10, ncol=4)
diagSigma <- matrix(0, ncol=nj, nrow=nj)
diag(diagSigma) <- NA
mydata <- list(
  nj = nj,
  m = m,
  id = matrix(1:nrow(datlong), ncol=nj, nrow=m, byrow=T),
  u=u,
  Nt = 2,
  Y = as.matrix(datlong[,3]),
  diagSigma = diagSigma
)

bcsm.model <- function(){

  for(j in 1:m){
    # data model
    Y[id[j,],] ~ dmnorm.vcov(beta[], omega[,])
  }
  # factor loadings
  for(i in 1:nj){
    lambda[i] ~ dnorm(0, 1)
  }
  # set u[,2] to factor loadings
  u[1:nj,2] <- lambda[]

  # item intercepts
  for(i in 1:nj){
    beta[i] ~ dnorm(0, 0.1)
  }

  # error variance
  for(i in 1:nj){
    sigma[i] ~ dunif(0,3)
    diagSigma[i,i] <- sigma[i]
  }
  SigmaInv[1, 1:nj, 1:nj] <- inverse(diagSigma[1:nj, 1:nj] )

  # =========================
  # layers
  # check value of theta for approx 0 value
  for(t in 2:(Nt+1)){
    LAM[t] <- ifelse(abs(theta[t-1]) < 1e-5, 1e5, 1/theta[t-1])
    d[t] <- LAM[t] + t(u[,t])%*%SigmaInv[t-1,1:nj, 1:nj]%*%u[,t]
    # update layer (t) inverse
    A[1:nj, t] <- SigmaInv[t-1,1:nj, 1:nj]%*%u[1:nj,t]
    for(i in 1:nj){
      for(j in 1:nj){
        B[t,i,j] <- A[i,t]*u[j,t]
      }
    }

    C[t,1:nj,1:nj] <- B[t, , ]%*%SigmaInv[t-1, , ]
    SigmaInv[t, 1:nj, 1:nj] <- SigmaInv[t-1, 1:nj, 1:nj] - (C[t,,])/d[t]
    # update truncation parameter
    tr[t-1] <- - 1/(t(u[,t])%*%SigmaInv[t-1, 1:nj, 1:nj]%*%u[,t])
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
  # update full covariance matrix
  omega[1:nj , 1:nj] <- diagSigma + utsum

}


# vector of all parameters to save
# exclude fixed lambda since it throws an error in
# in the GRB plot
param_save <- c("lambda", "sigma","beta","theta", "tr")

# helps to start all factor loadings as 1s
# makes the first iteration estimable
start_values <- list(
  list("lambda"=rep(1, nj)),
  list("lambda"=rep(1, nj))
)

# fit model
fit0 <- jags(
  model.file=bcsm.model,
  data=mydata,
  inits=start_values,
  parameters.to.save = param_save,
  n.iter=1000,n.chains = 2)

print(fit0)
fit <- fit0

fit <- R2jags::autojags(fit0,n.iter = 50000, n.thin = 20, Rhat = 1.01)
print(fit)

jags.mcmc <- as.mcmc(fit)
#jags.mcmc <- as.mcmc(fit$BUGSoutput$sims.matrix)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
plot.data <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(plot.data) <- c("chain", "iter", a)

mcmc_areas(plot.data,pars = "theta[1]", prob=0.8)
mcmc_areas(plot.data,pars = "theta[2]", prob=0.8)
mcmc_areas(plot.data,regex_pars = "beta", prob=0.8)
mcmc_areas(plot.data,regex_pars = "sigma", prob=0.8)

mcmc_acf(plot.data, regex_pars = "theta")
mcmc_acf(plot.data, regex_pars = "beta")
mcmc_acf(plot.data, regex_pars = "sigma")

mcmc_pairs(plot.data, regex_pars = "theta")

mcmc_trace(plot.data, regex_pars = "theta")
mcmc_trace(plot.data, regex_pars = "beta")
mcmc_trace(plot.data, regex_pars = "sigma")

gelman.plot(jags.mcmc)
