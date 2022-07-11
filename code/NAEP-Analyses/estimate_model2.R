# estimate_model2.R
# ===================================== #
# Math Motivation Study
# ===================================== #
# Created 
#   by: R. Noah Padgett
#   on: 2021-07-27
#
# Last Editted
#   by: R. Noah Padgett
#   on: 2021-08-10
# ===================================== #
# Purpose: 
# Estimate model 2 for results
#
# ===================================== #
# Load data
source("scripts/get_data.R")
# utility functions
source("scripts/get_utility_functions.R")


# Model 2: IFA with RT 
#This model is similar to the BL-IRT model for jointly modeling item responses and response times (Molenaar et al., 2015/2021).

mydata <- list(
  y   = as.matrix(sdat_full[,mathIdent]),
  lrt = pseudo_log(as.matrix(sdat_full[,colnames(sdat_full) %like% "RT_"])),
  N   = nrow(sdat_full),
  nit = length(mathIdent)
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
      # log-RT model
      dev[n,i]<-lambda[i]*(eta[n] - (tau[i,1]+tau[i,2]+tau[i,3]+tau[i,4])/2)
      mu.rt[n,i]<- icept[i] - speed[n] - rho * abs(dev[n,i])
      lrt[n,i]~dnorm(mu.rt[n,i],prec[i])
    }
  }
  ### Priors
  # person parameters
  for(n in 1:N){
    eta[n] ~ dnorm(0, 1) # latent ability
    speed[n]~dnorm(sigma.ts*eta[n],prec.s)  # latent speed
  }
  sigma.ts ~ dnorm(0, 0.1)
  prec.s ~ dt(1, 2, 1);T(0,)
  for(i in 1:nit){
    # lrt parameters
    icept[i]~dnorm(0,.1)
    icept.exp[i] = exp(icept[i])
    prec[i]~dt(1, 2, 1);T(0,)
    # Thresholds
    tau[i, 1] ~ dnorm(0.0,0.1)
    tau[i, 2] ~ dnorm(0, 0.1);T(tau[i, 1],)
    tau[i, 3] ~ dnorm(0, 0.1);T(tau[i, 2],)
    tau[i, 4] ~ dnorm(0, 0.1);T(tau[i, 3],)
    # factor Loadings
    lambda[i] ~ dnorm(0,.1);T(0,)
    lambda.std[i] = lambda[i]/pow(theta[i], 0.5)
    # latent response variances
    theta[i] = 1 + pow(lambda[i],2)
  }
  #inverted-U effect for item i
  rho~dnorm(0,.5);I(0,)
  # important parameters
  sigma.t <- pow(prec.s, -1) + pow(sigma.ts, 2)
  cor.ts <- sigma.ts/(pow(sigma.t,0.5))
  #reliability
  num = lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]
  omega.reli = (pow(num,2))/( pow(num,2) + pow(5,2))
}# End Model

jags.params <- c(
  # ability measurement model
  "tau", "lambda","lambda.std", "phi", "eta", "theta","pi",
  # speed measurement parameters
  "rho", "icept","prec.s", "sigma.t", "sigma.ts", "cor.ts","icept.exp",
  # Person-Item distance & # expected response time 
    # Plotting Values
    "mu.rt", "dev", "omega.reli")

# fit model
fit2 <- fit <- jags(
  model.file=jags.model
  , data=mydata
  , parameters.to.save = jags.params
  , n.iter   = 10000
  , n.chains = 4
  , n.thin = 5
)
# save fitted model
save(fit2,file = paste0(outputdir,"model2_fit.RData"))

load(paste0(outputdir,"model2_fit.RData"))

fit2$BUGSoutput$summary[!(rownames(fit2$BUGSoutput$summary) %like% "eta" | 
                           rownames(fit2$BUGSoutput$summary) %like% "dev" | 
                           rownames(fit2$BUGSoutput$summary) %like% "mu.rt"| 
                           rownames(fit2$BUGSoutput$summary) %like% "pi"), ]
load("misc/model2_fit.RData")
write.table(fit2$BUGSoutput$summary, file="output/posterior_summaries_model2.txt", sep = "\t")

# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit2)
fit.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit.mcmc) <- c("chain", "iter", colnames(as.data.frame(jags.mcmc[[1]])))

# Posterior Summary
# Density
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "tau", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-dens-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,pars = paste0("lambda.std[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-dens-lambdastd.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "theta", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "rho", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-rho.pdf"), width=3.5, height=5, units="in")
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "icept", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-icept.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "prec.s", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-prec.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "sigma.t", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-sigmat.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "cor.ts", prob = 0.8); ggsave(paste0(outputdir,"model2-mcmc-dens-corts.pdf"), width=3.5, height=5, units="in")

# posterior autocorrelation
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "tau"); ggsave(paste0(outputdir,"model2-mcmc-acf-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-acf-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit.mcmc,pars = paste0("theta[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-acf-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "rho"); ggsave(paste0(outputdir,"model2-mcmc-acf-rho.pdf"), width=3, height=5, units="in")
bayesplot::mcmc_acf(fit.mcmc,pars = paste0("icept[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-acf-icept.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "sigma.t"); ggsave(paste0(outputdir,"model2-mcmc-acf-sigmat.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "cor.ts"); ggsave(paste0(outputdir,"model2-mcmc-acf-corts.pdf"), width=3, height=5, units="in")
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "prec.s"); ggsave(paste0(outputdir,"model2-mcmc-acf-prec.pdf"), width=6, height=10, units="in")

# posterior traceplots
bayesplot::mcmc_trace(fit.mcmc,pars = paste0("theta[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-trace-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-trace-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "tau"); ggsave(paste0(outputdir,"model2-mcmc-trace-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "rho"); ggsave(paste0(outputdir,"model2-mcmc-trace-rho.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,pars = paste0("icept[",1:5,"]")); ggsave(paste0(outputdir,"model2-mcmc-trace-icept.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "sigma.t"); ggsave(paste0(outputdir,"model2-mcmc-trace-sigmat.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "cor.ts"); ggsave(paste0(outputdir,"model2-mcmc-trace-corts.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "prec.s"); ggsave(paste0(outputdir,"model2-mcmc-trace-prec.pdf"), width=6, height=10, units="in")

# posterior GRB convergence
fit.mcmc.ggs <- ggmcmc::ggs(jags.mcmc)
ggmcmc::ggs_grb(fit.mcmc.ggs, family="theta"); ggsave(paste0(outputdir,"model2-mcmc-grb-theta.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="tau"); ggsave(paste0(outputdir,"model2-mcmc-grb-tau.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="lambda"); ggsave(paste0(outputdir,"model2-mcmc-grb-lambda.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="rho"); ggsave(paste0(outputdir,"model2-mcmc-grb-rho.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="icept"); ggsave(paste0(outputdir,"model2-mcmc-grb-icept.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="sigma.t"); ggsave(paste0(outputdir,"model2-mcmc-grb-sigmat.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="cor.ts"); ggsave(paste0(outputdir,"model2-mcmc-grb-corts.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="prec.s"); ggsave(paste0(outputdir,"model2-mcmc-grb-prec.pdf"), width=6, height=10, units="in")

# posterior predictive distribution
Niter <- 200
fit2$model$recompile()
fit2.extra <- rjags::jags.samples(fit2$model, variable.names = "pi", n.iter = Niter)
N <- fit2$model$data()[["N"]]
nit <- C <- 5
n <- i <- iter <- 1
y.prob.ppc <- array(dim=c(Niter*2, nit, C))
for(iter in 1:(Niter*2)){
  iter.i <- ifelse(iter > Niter, iter - Niter, iter)
  y <- matrix(nrow=N, ncol=nit)
  for(i in 1:nit){
    for(n in 1:N){
      chain <- ifelse(iter <= Niter, 1, 2)
      y[n,i] <- sample(1:5, 1, prob = fit2.extra$pi[n,i, 1:5, iter.i, chain])
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

ymat <- fit2$model$data()[["y"]]
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
ggsave(paste0(outputdir,"fig3-model2-ppc.png"),plot = p, width=7, height=5, units="in")
ggsave(paste0(outputdir,"fig3-model2-ppc.pdf"),plot = p, width=7, height=5, units="in")



# Evaluate relationship between Eta-item distance and RT
lrt <- pseudo_log(as.matrix(sdat_full[,colnames(sdat_full) %like% "RT_"]))
person <- colMeans(fit.mcmc)[paste0("eta[",1:N,"]")]
plot.dat <- data.frame(
  id = 1:N,
  person,
  item = paste0("item_",1),
  PIdiff = colMeans(fit.mcmc)[paste0("dev[",1:N,",1]")],
  mu.rt = colMeans(fit.mcmc)[paste0("mu.rt[",1:N,",1]")],
  lrt = lrt[,1]
)
for(i in 2:5){
  p1 <-  data.frame(
    id = 1:N,
    person,
    item = paste0("item_",i),
    PIdiff = colMeans(fit.mcmc)[paste0("dev[",1:N,",",i,"]")],
    mu.rt = colMeans(fit.mcmc)[paste0("mu.rt[",1:N,",",i,"]")],
    lrt = lrt[,i]
  )
  plot.dat <- full_join(plot.dat, p1)
}

cols = c("Parabolic fit"="red", "GAM fit"="blue")
p <- ggplot(plot.dat, aes(x=PIdiff, mu.rt, group=item))+
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
ggsave(plot=p, filename = "figs/model_2_PIdiff_RT-murt.pdf", width = 7, heigh=5, units="in")

p <- ggplot(plot.dat, aes(x=PIdiff, lrt, group=item))+
  geom_point(alpha=0.25)+
  geom_smooth(se=F,formula = y ~ poly(x,2), aes(color="Parabolic fit"))+
  geom_smooth(se=F,aes(color='GAM fit'))+
  facet_wrap(.~item, nrow = 1)+
  scale_color_manual(values=cols, name="Fitted Regression Line")+
  labs(x="Person-Item Distance", y="log Response Time")+
  theme_classic()+
  theme(
    legend.position = "bottom"
  )
p
ggsave(plot=p, filename = "figs/model_2_PIdiff_RT-lrt.pdf", width = 7, heigh=5, units="in")

