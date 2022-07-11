# estimate_model4.R
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
# Estimate model 4 for results
#
# ===================================== #
# Load data
source("scripts/get_data.R")
# utility functions
source("scripts/get_utility_functions.R")

# Model 4: IFA with RT only
mydata <- list(
  y   = as.matrix(sdat_full[,mathIdent]),
  lrt = as.matrix(sdat_full[,colnames(sdat_full) %like% "RT_"]),
  N   = nrow(sdat_full),
  nit = length(mathIdent),
  C   = 5
)
jags.model <- function(){
  ### Model
  for(n in 1:N){
    for(i in 1:nit){
      # data model
      y[n,i] ~ dcat(omega[n,i, ])
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
      # MISCLASSIFICATION MODEL
      for(c in 1:C){
        # generate informative prior for misclassificaiton
        #   parameters based on RT
        for(ct in 1:C){
          alpha[n,i,ct,c] <- ifelse(c == ct,
                                    ilogit(lrt[n,i])*10,
                                    (1/(C-1))*(1-ilogit(lrt[n,i]))*10
          )
        }
        # sample misclassification parameters using the informative priors
        gamma[n,i,c,1:C] ~ ddirch(alpha[n,i,c,1:C])
        # observed category prob (Pr(y=c))
        omega[n,i, c] = gamma[n,i,c,1]*pi[n,i,1] + 
          gamma[n,i,c,2]*pi[n,i,2] + 
          gamma[n,i,c,3]*pi[n,i,3] + 
          gamma[n,i,c,4]*pi[n,i,4] + 
          gamma[n,i,c,5]*pi[n,i,5]
      }
    }
  }
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
    lambda[i] ~ dnorm(0,.1);T(0,)
    lambda.std[i] = lambda[i]/pow(theta[i], 0.5)
    # latent response variances
    theta[i] = 1 + pow(lambda[i],2)
  }
  #reliability
  num = lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]
  omega.reli = (pow(num,2))/( pow(num,2) + pow(5,2))
}# End Model

# vector of all parameters to save
jags.params <- c(
  # ability measurement model
  "tau", "lambda","lambda.std", "eta", "theta", "omega.reli", "gamma")

# fit4 model
fit4 <- jags(
    model.file=jags.model
  , data=mydata
  , parameters.to.save = jags.params
  , n.iter   = 10000
  , n.chains = 4
  , n.thin = 5
)
#save fitted model
save(fit4,file = paste0(outputdir,"model4_fit.RData"))
load(paste0(outputdir,"model4_fit.RData"))
# print results
fit4$BUGSoutput$summary[
  !(rownames(fit4$BUGSoutput$summary) %like% "eta" |  
      rownames(fit4$BUGSoutput$summary) %like% "gamma"), ]

load("misc/model4_fit.RData")
write.table(fit4$BUGSoutput$summary, file="output/posterior_summaries_model4.txt", sep = "\t")


# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit4)
fit4.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit4.mcmc) <- c("chain", "iter", colnames(as.data.frame(jags.mcmc[[1]])))

# Posterior Summary
# Density
bayesplot::mcmc_areas(fit4.mcmc, regex_pars="tau", prob=0.8); ggsave(paste0(outputdir,"model4-mcmc-dens-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit4.mcmc,pars = paste0("lambda[",1:5,"]"), prob = 0.8); ggsave(paste0(outputdir,"model4-mcmc-dens-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit4.mcmc,regex_pars = "lambda.std", prob = 0.8); ggsave(paste0(outputdir,"model4-mcmc-dens-lambda-std.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit4.mcmc,regex_pars = "theta", prob = 0.8); ggsave(paste0(outputdir,"model4-mcmc-dens-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_areas(fit4.mcmc,pars = paste0("eta[",1:5,"]"), prob = 0.8); ggsave(paste0(outputdir,"model4-mcmc-dens-eta1-5.pdf"), width=6, height=10, units="in")
# posterior autocorrelation
bayesplot::mcmc_acf(fit4.mcmc,regex_pars = "tau"); ggsave(paste0(outputdir,"model4-mcmc-acf-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit4.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model4-mcmc-acf-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit4.mcmc,regex_pars = "theta"); ggsave(paste0(outputdir,"model4-mcmc-acf-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_acf(fit4.mcmc,pars = paste0("eta[",1:5,"]")); ggsave(paste0(outputdir,"model4-mcmc-acf-eta.pdf"), width=6, height=10, units="in")
# posterior traceplots
bayesplot::mcmc_trace(fit4.mcmc,regex_pars = "theta"); ggsave(paste0(outputdir,"model4-mcmc-trace-theta.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit4.mcmc,pars = paste0("lambda[",1:5,"]")); ggsave(paste0(outputdir,"model4-mcmc-trace-lambda.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit4.mcmc,regex_pars = "tau"); ggsave(paste0(outputdir,"model4-mcmc-trace-tau.pdf"), width=6, height=10, units="in")
bayesplot::mcmc_trace(fit4.mcmc,pars = paste0("eta[",1:5,"]")); ggsave(paste0(outputdir,"model4-mcmc-trace.pdf"), width=6, height=10, units="in")

# posterior GRB convergence
fit4.mcmc.ggs <- ggmcmc::ggs(jags.mcmc)
ggmcmc::ggs_grb(fit4.mcmc.ggs, family="theta"); ggsave(paste0(outputdir,"model4-mcmc-grb-theta.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit4.mcmc.ggs, family="tau"); ggsave(paste0(outputdir,"model4-mcmc-grb-tau.pdf"), width=6, height=10, units="in")
ggmcmc::ggs_grb(fit4.mcmc.ggs, family="lambda"); ggsave(paste0(outputdir,"model4-mcmc-grb-lambda.pdf"), width=6, height=10, units="in")


# Posterior Predictive Check
Niter <- 1000
fit4$model$recompile()
fit4.extra <- rjags::jags.samples(fit4$model, variable.names = "omega", n.iter = Niter, thin = 5)
Niter <- 200
N <- fit4$model$data()[["N"]]
nit <- C <- 5
n <- i <- iter <- 1
y.prob.ppc <- array(dim=c(Niter*4, nit, C))
for(iter in 1:(Niter*4)){
  iter.i <- iter
  if(iter > Niter){
    iter.i <- iter - Niter
  }
  if(iter > 2*Niter){
    iter.i <- iter - 2*Niter
  }
  if(iter > 3*Niter){
    iter.i <- iter - 3*Niter
  }
  
  y <- matrix(nrow=N, ncol=nit)
  for(i in 1:nit){
    for(n in 1:N){
      chain <- 1
      if(iter > Niter){
        chain <- 2
      }
      if(iter > 2*Niter){
        chain <- 3
      }
      if(iter > 3*Niter){
        chain <- 4
      }
      y[n,i] <- sample(1:5, 1, prob = fit4.extra$omega[n,i, 1:5, iter.i, chain])
    }
    y.prob.ppc[iter,i,1] <- sum(y[,i]==1)/N
    y.prob.ppc[iter,i,2] <- sum(y[,i]==2)/N
    y.prob.ppc[iter,i,3] <- sum(y[,i]==3)/N
    y.prob.ppc[iter,i,4] <- sum(y[,i]==4)/N
    y.prob.ppc[iter,i,5] <- sum(y[,i]==5)/N
  }
}

yppcmat <- matrix(c(y.prob.ppc), ncol=1)
z <- expand.grid(1:(Niter*4), 1:5, 1:5)
yppcmat <- data.frame(  iter = z[,1], nit=z[,2], C=z[,3], yppc = yppcmat)

ymat <- fit4$model$data()[["y"]]
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
ggsave(paste0(outputdir,"fig3-model4-ppc.png"),plot = p, width=7, height=5, units="in")
ggsave(paste0(outputdir,"fig3-model4-ppc.pdf"),plot = p, width=7, height=5, units="in")


# compute average misclassification matrix parameters
i <- c <- cc <- 1
gammaList <- list()
for(i in 1:nit){
  i<-5
  gamma <- matrix(nrow=5, ncol=5)
  for(c in 1:C){
    for(cc in 1:C){
      gamma[c,cc] <- mean(colMeans(fit4.mcmc)[paste0("gamma[",1:N,",",i,",",c,",",cc,"]")])
    }
  }
  gammaList[[i]] <- gamma
}
gammaList

# plot dist of gammas
plot.dat <- t(matrix(colMeans(fit4.mcmc)[colnames(fit4.mcmc) %like% "gamma"], nrow=125))
A <- expand.grid(1:5, 1:5, 1:5) %>% as.data.frame() %>%
  mutate(
    g = paste0("gamma_", Var1,"_",Var2,"_",Var3)
  )
colnames(plot.dat) <- A$g
plot.dat <- plot.dat %>%
  as.data.frame() %>%
  mutate(id=1:nrow(plot.dat)) %>%
  pivot_longer(
    cols=contains("gamma"),
    names_to = c("C", "Cp", "item"),
    names_prefix = "gamma_",
    names_sep = "_",
    values_to = "gamma"
  )

p <- plot.dat %>%
  ggplot(aes(x=gamma, color=item, fill=item)) + 
    geom_histogram(binwidth = 0.005) +
    facet_grid(C~Cp)+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    theme_classic()+
    theme(
      legend.position="bottom",
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
p  
ggsave(paste0(outputdir,"model4-gamma-dist.pdf"),plot = p, width=10, height=10, units="in")

library(ggdist)
p <- plot.dat %>%
  ggplot(aes(x=item,y=gamma, fill=item)) + 
    ggdist::stat_halfeye(
      adjust=2, justification=-0.1,.width=0, point_colour=NA,
      normalize="all", fill="grey75"
    ) +
    geom_boxplot(
      width=.15, outlier.color = NA, alpha=0.5
    )+
    scale_color_grey()+
    scale_fill_grey()+
    facet_grid(C~Cp)+
    labs(x=NULL,y="Misclassification Probability (gamma)")+
    theme_bw()+
    theme(
      legend.position="bottom",
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
p  
ggsave(paste0(outputdir,"model4-gamma-dist-v2.pdf"),plot = p, width=10, height=10, units="in")
# Evaluate relationship between Eta-item distance and RT
lrt<- fit4$model$data()[["lrt"]]
person <- colMeans(fit4.mcmc)[paste0("eta[",1:N,"]")]
plot.dat <- data.frame(
  id = 1:N,
  person,
  item = paste0("item_",1),
  PIdiff = person - mean(colMeans(fit4.mcmc)[paste0("tau[1,",1:4,"]")]),
  lrt = lrt[,1]
)
for(i in 2:5){
  p1 <-  data.frame(
    id = 1:N,
    person,
    item = paste0("item_",i),
    PIdiff = person - mean(colMeans(fit4.mcmc)[paste0("tau[",i,",",1:4,"]")]),
    lrt = lrt[,i]
  )
  plot.dat <- full_join(plot.dat, p1)
}

cols = c("Parabolic fit"="red", "GAM fit"="blue")
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
ggsave(plot=p, filename = "output/model_4_PIdiff_RT.pdf", width = 7, heigh=5, units="in")

