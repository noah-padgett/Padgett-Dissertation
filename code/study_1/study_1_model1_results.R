# start of study 1 sim

# ===================================== #
# study_1_generate_data.R
# ===================================== #
# Padgett - Dissertation
# Created
#   on: 2022-01-06
#   by: R. Noah Padgett
# Last Editted
#   on: 2022-01-06
#   by: R. Noah Padgett
# ===================================== #
# Purpose: Generate data for study 1
# ===================================== #

# Load packages & utility functions
source("code/load_packages.R")
source("code/load_utility_functions.R")
# environment options
options(scipen = 999, digits=3)
# generate data for study 1
source("code/study_1/study_1_generate_data.R")


# summarise observed data

d1 <- sim.data %>%
  as.data.frame() %>%
  select(contains("y"))%>%
  mutate(id = 1:n()) %>%
  pivot_longer(
    cols=contains("y"),
    names_to = c("item"),
    values_to = "Response"
  ) %>%
  mutate(
    item = ifelse(nchar(item)==4,substr(item, 3,4),substr(item, 3,3))
  )
d2 <- sim.data %>%
  as.data.frame() %>%
  select(contains("logt"))%>%
  mutate(id = 1:n()) %>%
  pivot_longer(
    cols=contains("logt"),
    names_to = c("item"),
    values_to = "Time"
  ) %>%
  mutate(
    item = ifelse(nchar(item)==7,substr(item, 6,7),substr(item, 6,6))
  )
dat <- left_join(d1, d2)

dat_sum <- dat %>%
  select(item, Response, Time) %>%
  group_by(item) %>%
  summarize(
    p1 = table(Response)[1]/n(),
    p2 = table(Response)[2]/n(),
    p3 = table(Response)[3]/n(),
    M1 = mean(Response, na.rm=T),
    Mt = mean(Time, na.rm=T),
    SDt = sd(Time, na.rm=T)
  )

colnames(dat_sum) <- c("Item", "Prop. R == 1", "Prop. R == 2", "Prop. R == 3", "Mean Response", "Mean Response Time", "SD Response Time")
dat_sum$Item <- paste0("item_", 1:10)

kable(dat_sum, format="html", digits=3) %>%
  kable_styling(full_width = T)

  # Save parameters
  jags.params <- c(
    "tau", "lambda", "theta", "reli.omega"
  )

  mydata <- list(
    y = sim.data$y,
    N = nrow(sim.data$y),
    nit=ncol(sim.data$y)
  )

  # Run model
  # Model 1

  model.fit <-  R2jags::jags(
    model=paste0(w.d,"/code/study_1/model_1.txt")
    , parameters.to.save = jags.params
    , data=mydata
    , n.chains = 4
    , n.burnin = 5000
    , n.iter = 10000
    , n.thin = 5
  )

  print(model.fit)

# extract for plotting
jags.mcmc <- as.mcmc(model.fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
fit.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit.mcmc) <- c("chain", "iter", a)
fit.mcmc.ggs <- ggmcmc::ggs(jags.mcmc) # for GRB plot

# Posterior Summary
# tau
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "tau", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "tau")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "tau")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="tau")

bayesplot::mcmc_areas(fit.mcmc,regex_pars = "lambda", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "lambda")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "lambda")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="lambda")

bayesplot::mcmc_areas(fit.mcmc,regex_pars = "theta", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "theta")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "theta")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="theta")

bayesplot::mcmc_areas(fit.mcmc,regex_pars = "reli.omega", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "reli.omega")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "reli.omega")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="reli.omega")


# extract omega posterior for results comparison
extracted_omega <- data.frame(model_1=fit.mcmc$reli.omega)



# Model 2
# Save parameters
  jags.params <- c(
    "tau", "lambda", "theta",
    "icept", "prec", "prec.s", "sigma.ts", "rho",
    "reli.omega"
  )

  mydata <- list(
    y = sim.data$y,
    lrt=sim.data$logt,
    N = nrow(sim.data$y),
    nit=ncol(sim.data$y)
  )

  # Run model
  # Model 2

  model.fit <-  R2jags::jags(
    model=paste0(w.d,"/code/study_1/model_2.txt")
    , parameters.to.save = jags.params
    , data=mydata
    , n.chains = 4
    , n.burnin = 5000
    , n.iter = 10000
    , n.thin = 5
  )

  print(model.fit)

# extract for plotting
jags.mcmc <- as.mcmc(model.fit)
a <- colnames(as.data.frame(jags.mcmc[[1]]))
fit.mcmc <- data.frame(as.matrix(jags.mcmc, chains=T, iters = T))
colnames(fit.mcmc) <- c("chain", "iter", a)
fit.mcmc.ggs <- ggmcmc::ggs(jags.mcmc) # for GRB plot

# Posterior Summary
# tau
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "tau", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "tau")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "tau")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="tau")
# lambda
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "lambda", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "lambda")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "lambda")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="lambda")
# theta
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "theta", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "theta")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "theta")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="theta")
# icept
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "icept", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "icept")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "icept")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="icept")
# prec
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "prec", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "prec")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "prec")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="prec")
# prec.s
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "prec.s", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "prec.s")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "prec.s")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="prec.s")
# sigma.ts
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "sigma.ts", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "sigma.ts")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "sigma.ts")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="sigma.ts")
# rho
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "rho", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "rho")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "rho")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="rho")

# reliability
bayesplot::mcmc_areas(fit.mcmc,regex_pars = "reli.omega", prob = 0.8)
bayesplot::mcmc_acf(fit.mcmc,regex_pars = "reli.omega")
bayesplot::mcmc_trace(fit.mcmc,regex_pars = "reli.omega")
ggmcmc::ggs_grb(fit.mcmc.ggs, family="reli.omega")

# extract omega posterior for results comparison
extracted_omega$model_2 <- fit.mcmc$reli.omega


