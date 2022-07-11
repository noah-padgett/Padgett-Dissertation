# estimate_model7.R
# ===================================== #
# Math Motivation Study
# ===================================== #
# Created 
#   by: R. Noah Padgett
#   on: 2021-07-27
#
# Last Editted
#   by: R. Noah Padgett
#   on: 2021-08-28
# ===================================== #
# Purpose: 
# Estimate model 7 for results
#
# ===================================== #
# Load data
source("scripts/get_data.R")
# utility functions
source("scripts/get_utility_functions.R")
# Model 7: Full IFA with Misclassification
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
          gamma[n,i,ct,c] <- ifelse(c == ct,
                                    ilogit(lrt[n,i]),
                                    (1/(C-1))*(1-ilogit(lrt[n,i]))
          )
        }
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
    eta[n] ~ dnorm(0, 1) # latent ability
  }

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
} # End Model

jags.params <- c(
  # ability measurement model
  "tau", "lambda","lambda.std",  "theta", "eta",
  # misclassification parameters
  "gamma"
  # misc. 
  , "omega.reli")

# fit7 model
fit7 <- jags(
  model.file=jags.model
  , data=mydata
  , parameters.to.save = jags.params
  , n.iter   = 10000
  , n.chains = 4
  , n.thin = 5
)
# save fitted model
save(fit7,file = "misc/model7_fit.RData")

#load(paste0(outputdir,"model6_fit.RData"))
# print results
fit7$BUGSoutput$summary[
  !(rownames(fit7$BUGSoutput$summary) %like% "eta" |  
      rownames(fit7$BUGSoutput$summary) %like% "gamma"), ]
write.table(fit7$BUGSoutput$summary, file="output/posterior_summaries_model7.txt", sep = "\t")
# extract posteriors for all chains
jags.mcmc <- as.mcmc(fit7)
fit7.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit7.mcmc) <- c("chain", "iter", colnames(as.data.frame(jags.mcmc[[1]])))

# Posterior Summary
# Density
bayesplot::mcmc_areas(fit7.mcmc,regex_pars = "tau", prob = 0.8)
bayesplot::mcmc_areas(fit7.mcmc,regex_pars = "lambda", prob = 0.8)
bayesplot::mcmc_areas(fit7.mcmc,regex_pars = "theta", prob = 0.8)

# posterior autocorrelation
bayesplot::mcmc_acf(fit7.mcmc,regex_pars = "tau")
bayesplot::mcmc_acf(fit7.mcmc,regex_pars = "lambda")
bayesplot::mcmc_acf(fit7.mcmc,regex_pars = "theta")

# posterior traceplots
bayesplot::mcmc_trace(fit7.mcmc,regex_pars = "theta")
bayesplot::mcmc_trace(fit7.mcmc,regex_pars = "lambda")
bayesplot::mcmc_trace(fit7.mcmc,regex_pars = "tau")

# posterior GRB convergence
fit7.mcmc.ggs <- ggmcmc::ggs(jags.mcmc)
ggmcmc::ggs_grb(fit7.mcmc.ggs, family="theta")
ggmcmc::ggs_grb(fit7.mcmc.ggs, family="tau")
ggmcmc::ggs_grb(fit7.mcmc.ggs, family="lambda")


# Posterior Predictive Check
Niter <- 200
fit7$model$recompile()
fit7.extra <- rjags::jags.samples(fit7$model, variable.names = "omega", n.iter = Niter)
N <- fit7$model$data()[["N"]]
nit <- C <- 5
n <- i <- iter <- 1
y.prob.ppc <- array(dim=c(Niter*2, nit, C))
for(iter in 1:(Niter*2)){
  iter.i <- ifelse(iter > Niter, iter - Niter, iter)
  y <- matrix(nrow=N, ncol=nit)
  for(i in 1:nit){
    for(n in 1:N){
      chain <- ifelse(iter <= Niter, 1, 2)
      y[n,i] <- sample(1:5, 1, prob = fit7.extra$omega[n,i, 1:5, iter.i, chain])
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

ymat <- fit7$model$data()[["y"]]
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
  facet_wrap(.~nit, nrow=1)+
  theme_classic()
p
ggsave(paste0(outputdir,"fig3-model7-ppc.png"),plot = p, width=7, height=5, units="in")
ggsave(paste0(outputdir,"fig3-model7-ppc.pdf"),plot = p, width=7, height=5, units="in")


# compute average misclassification matrix parameters
i <- c <- cc <- 1
gammaList <- list()
for(i in 1:nit){
  gamma <- matrix(nrow=5, ncol=5)
  for(c in 1:C){
    for(cc in 1:C){
      gamma[c,cc] <- mean(colMeans(fit7.mcmc)[paste0("gamma[",1:N,",",i,",",c,",",cc,"]")])
    }
  }
  gammaList[[i]] <- gamma
}
gammaList

plot.dat <- t(matrix(colMeans(fit7.mcmc)[colnames(fit7.mcmc) %like% "gamma"], nrow=125))
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
ggsave(paste0(outputdir,"model7-gamma-dist.pdf"),plot = p, width=10, height=10, units="in")

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
ggsave(paste0(outputdir,"model7-gamma-dist-v2.pdf"),plot = p, width=10, height=10, units="in")
# Evaluate relationship between Eta-item distance and RT
lrt <- pseudo_log(as.matrix(sdat_full[,colnames(sdat_full) %like% "RT_"]))
person <- colMeans(fit7.mcmc)[paste0("eta[",1:N,"]")]
plot.dat <- data.frame(
  id = 1:N,
  person,
  item = paste0("item_",1),
  PIdiff = person - mean(colMeans(fit7.mcmc)[paste0("tau[1,",1:4,"]")]),
  lrt = lrt[,1]
)
for(i in 2:5){
  p1 <-  data.frame(
    id = 1:N,
    person,
    item = paste0("item_",i),
    PIdiff = person - mean(colMeans(fit7.mcmc)[paste0("tau[",i,",",1:4,"]")]),
    lrt = lrt[,i]
  )
  plot.dat <- full_join(plot.dat, p1)
}

cols = c("Parabolic fit"="red", "GAM fit"="blue")
p <- ggplot(plot.dat, aes(x=PIdiff, lrt, group=item))+
  geom_point()+
  geom_smooth(se=F,formula = y ~ poly(x,2), aes(color="Parabolic fit"))+
  geom_smooth(se=F,aes(color='GAM fit'))+
  facet_wrap(.~item, scales = "free_x",nrow = 1)+
  scale_color_manual(values=cols, name="Fitted Regression Line")+
  labs(x="Person-Item Distance", y="log Response Time (model expected)")+
  theme_classic()+
  theme(
    legend.position = "bottom"
  )
p
ggsave(plot=p, filename = "figs/model_7_PIdiff_RT.pdf", width = 7, heigh=5, units="in")


