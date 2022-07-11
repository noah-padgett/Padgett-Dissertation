# estimate_model1.R
# ===================================== #
# Math Motivation Study
# ===================================== #
# Created 
#   by: R. Noah Padgett
#   on: 2021-07-27
#
# Last Editted
#   by: R. Noah Padgett
#   on: 2021-08-01
# ===================================== #
# Purpose: 
# Estimate model 1 for getting person
#   scores (eta) and factor model
#   parameters
#
# ===================================== #
# Load data
source("scripts/get_data.R")
# utility functions
source("scripts/get_utility_functions.R")

# =================================== #
# Model 1: Traditional IFA

# parameterization => conditional standardized factor
# Var(factor) = 1
# Var(Item Residual) = 1
# => Total Variance Item Latent Response Variable 
#   Var(Item Latent Response) = 1 + lambda^2

# set up estimation
mydata <- list(
  y = as.matrix(sdat_full[,mathIdent]),
  N = nrow(sdat_full),
  nit=length(mathIdent)
)
jags.model <- function(){
  ### Model
  for(n in 1:N){
    for(i in 1:nit){
      # data model
      y[n,i] ~ dcat(pi[n,i, ])
      
      # LRV
      ystar[n,i] ~ dnorm(lambda[i]*eta[n], 1)
      
      # Pr(nu = 5)
      pi[n,i,5] = phi(ystar[n,i] - tau[i,4])
      # Pr(nu = 4)
      pi[n,i,4] = phi(ystar[n,i] - tau[i,3]) - phi(ystar[n,i] - tau[i,4])
      # Pr(nu = 3)
      pi[n,i,3] = phi(ystar[n,i] - tau[i,2]) - phi(ystar[n,i] - tau[i,3])
      # Pr(nu = 2)
      pi[n,i,2] = phi(ystar[n,i] - tau[i,1]) - phi(ystar[n,i] - tau[i,2])
      # Pr(nu = 1)
      pi[n,i,1] = 1 - phi(ystar[n,i] - tau[i,1])
    }
  }
  ### Priors
  # person parameters
  for(n in 1:N){
    eta[n] ~ dnorm(0, 1)
  }
  # item parameters
  for(i in 1:nit){
    # Thresholds
    tau[i, 1] ~ dnorm(0.0,0.1)
    tau[i, 2] ~ dnorm(0, 0.1);T(tau[i, 1],)
    tau[i, 3] ~ dnorm(0, 0.1);T(tau[i, 2],)
    tau[i, 4] ~ dnorm(0, 0.1);T(tau[i, 3],)
    # factor Loadings
    lambda[i] ~ dgamma(2,2);T(0,)
    lambda.std[i] = lambda[i]/pow(theta[i], 0.5)
    # latent response variances
    theta[i] = 1 + pow(lambda[i],2)
  }
  #reliability
  num = lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]
  omega.reli = (pow(num,2))/( pow(num,2) + pow(5,2))
}# End Model

# vector of all parameters to save
param_save <- c("lambda","lambda.std", "tau","theta", "eta", "omega.reli")

# fit1 model
fit1 <- jags(
  model.file=jags.model
  , data=mydata
  , parameters.to.save = jags.params
  , n.iter   = 10000
  , n.chains = 4
  , n.thin = 5
)

load("misc/model1_fit.RData")

print(fit1)
fit1$BUGSoutput$summary[!(rownames(fit1$BUGSoutput$summary) %like% "eta" | 
                            rownames(fit1$BUGSoutput$summary) %like% "pi"), ]
load("misc/model1_fit.RData")
write.table(fit1$BUGSoutput$summary, file="output/posterior_summaries_model1.txt", sep = "\t")
# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit1)
fit1.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit1.mcmc) <- c("chain", "iter", colnames(as.data.frame(jags.mcmc[[1]])))

# Posterior Summary
# Density
bayesplot::mcmc_areas(fit1.mcmc, regex_pars="tau", prob=0.8); ggsave(paste0(outputdir,"model1-mcmc-dens-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit1.mcmc,pars = paste0("lambda[",1:5,"]"), prob = 0.8); ggsave(paste0(outputdir,"model1-mcmc-dens-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit1.mcmc,regex_pars = "lambda.std", prob = 0.8); ggsave(paste0(outputdir,"model1-mcmc-dens-lambda-std.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit1.mcmc,regex_pars = "theta", prob = 0.8); ggsave(paste0(outputdir,"model1-mcmc-dens-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit1.mcmc,pars = paste0("eta[",1:5,"]"), prob = 0.8); ggsave(paste0(outputdir,"model1-mcmc-dens-eta1-5.pdf"), width=6, height=10, units="in")
# posterior autocorrelation
bayesplot::mcmc_acf(fit1.mcmc,regex_pars = "tau"); ggsave(paste0(outputdir,"model1-mcmc-acf-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit1.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model1-mcmc-acf-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit1.mcmc,regex_pars = "theta"); ggsave(paste0(outputdir,"model1-mcmc-acf-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit1.mcmc,pars = paste0("eta[",1:5,"]")); ggsave(paste0(outputdir,"model1-mcmc-acf-eta.pdf"), width=6, height=10, units="in")
# posterior traceplots
bayesplot::mcmc_trace(fit1.mcmc,regex_pars = "theta"); ggsave(paste0(outputdir,"model1-mcmc-trace-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit1.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model1-mcmc-trace-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit1.mcmc,regex_pars = "tau"); ggsave(paste0(outputdir,"model1-mcmc-trace-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit1.mcmc,pars = paste0("eta[",1:5,"]")); ggsave(paste0(outputdir,"model1-mcmc-trace.pdf"), width=6, height=10, units="in")

# posterior GRB convergence
fit1.mcmc.ggs <- ggmcmc::ggs(jags.mcmc)
ggmcmc::ggs_grb(fit1.mcmc.ggs, family="theta"); ggsave(paste0(outputdir,"model1-mcmc-grb-theta.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit1.mcmc.ggs, family="tau"); ggsave(paste0(outputdir,"model1-mcmc-grb-tau.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit1.mcmc.ggs, family="lambda"); ggsave(paste0(outputdir,"model1-mcmc-grb-lambda.pdf"), width=6, height=10, units="in")

# Posterior Predictive Check
Niter <- 200
fit1$model$recompile()
fit1.extra <- rjags::jags.samples(fit1$model, variable.names = "pi", n.iter = Niter)
N <- fit1$model$data()[[1]]
nit <- C <- 5
n <- i <- iter <- 1
y.prob.ppc <- array(dim=c(Niter*2, nit, C))
for(iter in 1:(Niter*2)){
  iter.i <- ifelse(iter > Niter, iter - Niter, iter)
  y <- matrix(nrow=N, ncol=nit)
  for(i in 1:nit){
    for(n in 1:N){
      chain <- ifelse(iter <= Niter, 1, 2)
      y[n,i] <- sample(1:5, 1, prob = fit1.extra$pi[n,i, 1:5, iter.i, chain])
    }
    y.prob.ppc[iter,i,1] <- sum(y[,i]==1)/N
    y.prob.ppc[iter,i,2] <- sum(y[,i]==2)/N
    y.prob.ppc[iter,i,3] <- sum(y[,i]==3)/N
    y.prob.ppc[iter,i,4] <- sum(y[,i]==4)/N
    y.prob.ppc[iter,i,5] <- sum(y[,i]==5)/N
  }
}

yppcmat <- matrix(c(y.prob.ppc), ncol=1)
z <- expand.grid(1:(Niter*2), 1:5, 1:5)
yppcmat <- data.frame(  iter = z[,1], nit=z[,2], C=z[,3], yppc = yppcmat)

ymat <- fit1$model$data()[[3]]
y.prob <- matrix(ncol=C, nrow=nit)
for(i in 1:nit){
  for(c in 1:C){
    y.prob[i,c] <- sum(ymat[,i]==c)/N
  }
}
yobsmat <- matrix(c(y.prob), ncol=1)
z <- expand.grid(1:5, 1:5)
yobsmat <- data.frame(nit=z[,1], C=z[,2], yobs = yobsmat)
plot.ppc <- full_join(yppcmat, yobsmat)

p <- plot.ppc %>%
  mutate(C    = as.factor(C),
         item = nit) %>%
  ggplot()+
  geom_boxplot(aes(x=C,y=y.prob.ppc), outlier.colour = NA)+
  geom_point(aes(x=C,y=yobs), color="red")+
  scale_color_grey()+
  facet_wrap(.~nit, nrow=1)+
  theme_classic()
p
ggsave(paste0(outputdir,"fig3-model1-ppc.png"),plot = p, width=7, height=5, units="in")
ggsave(paste0(outputdir,"fig3-model1-ppc.pdf"),plot = p, width=7, height=5, units="in")

# Evaluate relationship between Eta and RT
lrt <- pseudo_log(as.matrix(sdat_full[,colnames(sdat_full) %like% "RT_"]))
person <- colMeans(fit1.mcmc)[paste0("eta[",1:N,"]")]

plot.dat <- data.frame(
  id = 1:N,
  person,
  item = paste0("item_",1),
  PIdiff = person - mean(colMeans(fit1.mcmc)[paste0("tau[1,",1:4,"]")]),
  lrt = lrt[,1]
)
for(i in 2:5){
  p1 <-  data.frame(
    id = 1:N,
    person,
    item = paste0("item_",i),
    PIdiff = person - mean(colMeans(fit1.mcmc)[paste0("tau[",i,",",1:4,"]")]),
    lrt = lrt[,i]
  )
  plot.dat <- full_join(plot.dat, p1)
}

cols = c("Parabolic fit"="red", "GAM fit"="blue")
p <- ggplot(plot.dat, aes(x=PIdiff, lrt, group=item))+
  geom_point(alpha=0.25)+
  geom_smooth(se=F,formula = y ~ poly(x,2), aes(color="Parabolic fit"))+
  geom_smooth(se=F,aes(color='GAM fit'))+
  facet_wrap(.~item,nrow = 1)+
  scale_color_manual(values=cols, name="Fitted Regression Line")+
  labs(x="Person-Item Distance", y="log Response Time")+
  theme_classic()+
  theme(
    legend.position = "bottom"
  )
p
ggsave(plot=p, filename = "output/model_1_PIdiff_RT.pdf", width = 7, heigh=5, units="in")




