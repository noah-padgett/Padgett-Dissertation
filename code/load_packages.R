# ============================================= #
# script: load_packages.R
# Project: PhD Dissertation Code
# Author(s): R.N. Padgett et al.
# ============================================= #
# Data Created: 2020-11-28
# Date Modified: 2021-04-14
# By: R. Noah Padgett
# ============================================= #
# Purpose:
# This R script is for loading all necessary
#   R packages
#
# No output - just loading packages into the
#   environment
# ============================================= #
# Set up directory and libraries
rm(list=ls())
# list of packages
packages <- c(
  "tidyverse", "readr", "patchwork",
  "tidyr","data.table", "dplyr","ggplot2",
  "coda", "ggmcmc", "bayesplot", "polycor", "lavaan",
  "kableExtra", "xtable",
  "diffIRT","eRm", "R2jags", "sirt", "lme4", 'runjags',
  "LaplacesDemon", "mvtnorm", "car"
)
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load packages
lapply(packages, library, character.only = TRUE)

w.d <- getwd()
