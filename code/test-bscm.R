# covariances structure modeling
source("code/load_packages.R")
# sample sizes
nj <- 4 # i
m <- 10 # j
# matrices
I <- diag(nrow=nj, ncol=nj)
J <- matrix(1, ncol=nj, nrow=nj)
# (co)variance parameters
d <- 1
sigma2 <- 1*d
tau.mu <- -sigma2/nj + 0.1
tau.sig <- 0.05*d
# function to generate RE from a truncated normal distr.
rnormt <- function(n, mu, s = 1, ll=-Inf, ul=Inf) {
  # range is a vector of two values
  F.a <- pnorm(ll, mean = mu, sd = s)
  F.b <- pnorm(ul, mean = mu, sd = s)
  # random uniform within truncated range
  u <- runif(n, min = F.a, max = F.b)
  # use the quantile "u" to obtain random draw
  qnorm(u, mean = mu, sd = s)
}
tau.re <- rnormt(m, tau.mu, tau.sig, -sigma2/nj, sigma2/nj) - tau.mu

# matrix covariance matrix prior to covariate adjustment
# expected value no random component
omega <- I*sigma2 + J*tau.mu
omega

# fixed effects
beta <- c(0.3, 0.2)

# data vector
dat <- matrix(nrow=nj*m, ncol=3)
i <- 1; j <- 1
for(j in 1:m){
  dat[i:(i+nj-1), 1] <- j
  dat[i:(i+nj-1), 2] <- X<- rnorm(nj, 0, 0.25) # X
  # compute error
  omegaj <- omega + J*tau.re[j]
  Ej <- mvtnorm::rmvnorm(1, sigma=omegaj)
  dat[i:(i+nj-1), 3] <- beta[1] + X*beta[2] + Ej# Y
  i <- i + nj
}

dat <- as.data.frame(dat)
colnames(dat) <- c("g","X", "Y")
dat$g <- as.factor(dat$g)

# plot
ggplot(dat, aes(x=X, y=Y, color=g))+
  geom_point()+
  geom_smooth(method="lm", se=F)+
  facet_wrap(g~.)

# lme4
fit <- lmer(Y ~ 1+ (1| g), data=dat)
summary(fit)

fit <- lmer(Y ~ 1+ X + (1+X| g), data=dat)
summary(fit)

for(j in 1:m){
  print(lm(Y ~ 1 + X, data=filter(dat, g==j)))
}

# jags

# model code
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

dat

mydata <- list(
  nj = nj,
  m = m,
  id = matrix(1:nrow(dat), ncol=nj, nrow=m, byrow=T),
  Imat = diag(nrow=nj, ncol=nj),
  Ivec = matrix(1, nrow=nj, ncol=1),
  Jmat = matrix(1, ncol=nj, nrow=nj),
  Y = dat[,3,drop=F] #, X = dat[,2,drop=F]
)


# initial values
start_values <- list(
  list("sigma2"=1,
       "tau.mu"=0,
       "tau.sig"=0.05,
       "tau"=rep(0,m),
       "beta"=c(0))
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
  n.iter=10000)

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

plot.data$log.tau.sig <- log(plot.data$tau.sig)
mcmc_pairs(plot.data, pars=c("tau.mu", "tau.sig", "sigma2"))
mcmc_pairs(plot.data, pars=c("tau.mu", "log.tau.sig", "sigma2"))

mcmc_acf(plot.data, pars=c("tau.mu", "tau.sig", "sigma2"))
mcmc_acf(plot.data, pars=paste0("tau[",1:10,"]"))


# prios
sim_prior <- function(nj=10){
  #inv.sigma <- runif(1, 0,3) #rgamma(1, 5, 10) # error precision
  sigma2 <- runif(1,0,3)
  lb <- -sigma2/nj + 0.001
  tau.mu <- -1
  i <- 1
  while(tau.mu < lb){
    tau.mu <- rnorm(1, 0, 0.5)
    i <- i + 1
    if(i == 100) break
  }
  tau.sig <- runif(1,0,0.2) # rgamma(1, 1, 20)#rlnorm(1, -4,1/(1.75**2))
  c(sigma2, tau.mu, tau.sig)
}
rprior <- function(n){
  sims <- matrix(nrow=n, ncol=3)
  sims <- t(replicate(n, sim_prior()))
  sims
}

prior.sims <- as.data.frame(rprior(10000))
colnames(prior.sims) <- c("sigma2", "tau.mu", "tau.sig")

cols <- c("Posterior"="#0072B2", "Prior"="#E69F00", "True"= "black")#"#56B4E9", "#E69F00" "#CC79A7"

# make plotting pieces
p1 <- ggplot()+
  geom_density(data=plot.data, aes(x=sigma, color="Posterior"))+
  geom_density(data=prior.sims, aes(x=sigma,color="Prior"),adjust=2)+
  geom_vline(aes(xintercept=1, color="True"))+
  scale_color_manual(values=cols, name=NULL)+
  lims(x=c(0, 3))+
  theme_classic()
p2 <- ggplot()+
  geom_density(data=plot.data, aes(x=tau.mu, color="Posterior"))+
  geom_density(data=prior.sims, aes(x=tau.mu,color="Prior"),adjust=2)+
  geom_vline(aes(xintercept=-0.075, color="True"))+
  scale_color_manual(values=cols, name=NULL)+
  lims(x=c(-0.2, 0.5))+
  theme_classic()
p3 <- ggplot()+
  geom_density(data=plot.data, aes(x=tau.sig, color="Posterior"))+
  geom_density(data=prior.sims, aes(x=tau.sig,color="Prior"),adjust=2)+
  geom_vline(aes(xintercept=0.01, color="True"))+
  scale_color_manual(values=cols, name=NULL)+
  lims(x=c(0, 0.5))+
  theme_classic()
p1 + p2 + p3 + plot_layout(guides="collect")


mcmc_pairs(prior.sims, pars=c("tau.mu", "tau.sig", "sigma"))
mcmc_pairs(plot.data, pars=c("tau.mu", "tau.sig", "sigma"))
