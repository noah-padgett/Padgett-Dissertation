# compare_omega.R
# ===================================== #
# Math Motivation Study
# ===================================== #
# Created
#   by: R. Noah Padgett
#   on: 2021-08-05
#
# Last Editted
#   by: R. Noah Padgett
#   on: 2021-08-28
# ===================================== #
# Purpose:
# Compare estimates of reliability
#
# ===================================== #
# Load data
source("scripts/get_data.R")
# utility functions
source("scripts/get_utility_functions.R")

# load fitted models
load("misc/model1_fit.RData")
load("misc/model2_fit.RData")
load("misc/model3_fit.RData")
load("misc/model4_fit.RData")
load("misc/model5_1_fit.RData")
load("misc/model6_fit.RData")
load("misc/model7_fit.RData")

# extract estimates
omega1 <- fit1$BUGSoutput$sims.matrix[, "omega.reli"]
omega2 <- fit2$BUGSoutput$sims.matrix[, "omega.reli"]
omega3 <- fit3$BUGSoutput$sims.matrix[, "omega.reli"]
omega4 <- fit4$BUGSoutput$sims.matrix[, "omega.reli"]
omega5 <- fit5$BUGSoutput$sims.matrix[, "omega.reli"]
omega6 <- fit6$BUGSoutput$sims.matrix[, "omega.reli"]
omega7 <- fit7$BUGSoutput$sims.matrix[, "omega.reli"]

rm(fit1, fit2, fit3, fit4, fit5, fit6, fit7)


dat <- data.frame(
  Model_1 = omega1,
  Model_2 = omega2,
  Model_3 = omega3,
  Model_4 = omega4,
  Model_5 = omega5,
  Model_6 = omega6,
  Model_7 = omega7
) %>%
  pivot_longer(
    cols=everything(),
    names_to = "model",
    values_to = "omega"
  )

describeBy(dat$omega,group = dat$model, mat=T)

p <- ggplot(dat, aes(x=model, y=omega))+
  geom_boxplot(outlier.alpha = 0)+
  lims(y=c(0.7, 1))+
  labs(x=NULL)+
  theme_classic()
p
ggsave(paste0(outputdir,"fig2-omega-compare-v1.2.png"),plot = p, width=5, height=3.5, units="in")
ggsave(paste0(outputdir,"fig2-omega-compare-v1.2.pdf"),plot = p, width=5, height=3.5, units="in")



library(ggdist)
p <- ggplot(dat, aes(x=model, y=omega))+
  ggdist::stat_halfeye(
    adjust=2, justification=-0.1,.width=0, point_colour=NA,
    normalize="all", fill="grey75"
  ) +
  geom_boxplot(
    width=.15, outlier.color = NA, alpha=0.5
  ) +
  labs(x=NULL)+
  lims(y=c(0.7, 1))+
  theme_classic()
p
ggsave(paste0(outputdir,"fig2-omega-compare-v2.2.png"),plot = p, width=5, height=3.5, units="in")
ggsave(paste0(outputdir,"fig2-omega-compare-v2.2.pdf"),plot = p, width=5, height=3.5, units="in")
