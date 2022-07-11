# ===================================== #
# cond_1.R
# ===================================== #
# Padgett - Dissertation
# Created
#   on: 2021-11-14
#   by: R. Noah Padgett
# Last Editted
#   on: 2021-12-06
#   by: R. Noah Padgett
# ===================================== #
# Purpose: Run condition 1 of MC Sim.
#   Condition 1:
#	N 	= 500
#	N Items = 5
#	N Cat   = 2
#	Nrep 	= 100
COND_NUM = 1
N_cat = 2
N_items = 5
N_persons = 500
N_rep = 100
# ===================================== #
# libraries
library(doParallel)
library(R2jags)
library(mvtnorm)



simulate_data_misclass <- function(paravec, tau=NULL){
  # NOTE: tau is a matrix[J, C-1] of item threshold parameters that possibly vary over items
  # useful functions
  invlogit <- function(x) {exp(x)/(1+exp(x))}
  logit <- function(x){log(x/(1-x))}
  # Generating Data
  N <- paravec[1] # number of respondents
  J <- paravec[2] # number of items
  C <- paravec[3] # number of response categories
  # ========================= #
  # latent person parameters
  etaCor <- paravec[4] # correlation between ability and speediness
  etasd <- paravec[5:6]
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
  lambda <- matrix(rep(paravec[7], J), ncol=1)
  # item latent residual variances
  theta <- c(1 - lambda**2)
  # item thresholds
  if(is.null(tau)){
    tau <- matrix(ncol=C-1, nrow=J)
    for(c in 1:(C-1)){
      if(c == 1){
        tau[,1] <- runif(J, -1, -0.33)
      }
      if(c > 1){
        tau[,c] <- tau[,c-1] + runif(J, 0.25, 1)
      }
    }
  }

  # latent item response
  ystar <- lambda%*%eta0
  ystar <- apply(ystar, 2, FUN = function(x){mvtnorm::rmvnorm(1, x, diag(theta, ncol=J, nrow=J))})
  # response time parameters (copied from Molenaar et al. 2021)
  nu <- matrix(rep(paravec[8], J), ncol=1)
  sigma.ei <- matrix(rep(paravec[9], J), ncol=1)
  rho1 <- paravec[10]
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
      # rescale to 0/1 if dichotomous items
      if(C == 2){
        Y[i,j] = Y[i,j]-1
      }
    }
  }
  # true_values <- list(eta0, eta1, lambda, nu, sigma.ei, tau, mulogt, ystar, theta, gamma, omega)
  # names(true_values) <- c("eta", "")
  sim_data <- list(Y, logt)
  names(sim_data) <- c("y", "logt")
  return(sim_data)

}


call_fx <- function(con)
{
  c<-con$c
  i<-con$i
  simulate_data_misclass <- con$fx_sim_dat
  # get data
  paravec <- c(
    N=con$N_persons, J = con$N_item, C = con$N_cat,
    etaCor = .23, etasd1 = 1, etasd2 = sqrt(0.1),
    lambda=0.9, nu=1.5, sigma.ei=0.25, rho1=0.1)
  sim.data <- simulate_data_misclass(paravec, tau=con$tau)

  # Save parameters (ignore misclass on kodaik)
  jags.params <- c(
    "tau", "lambda", "lambda.std", "theta",
    "icept", "prec", "prec.s", "sigma.ts", "rho",
    "reli.omega"
  )

  mydata <- list(
    y = sim.data$y,
    lrt = sim.data$logt,
    N = nrow(sim.data$y),
    nit=ncol(sim.data$y),
    ncat = paravec[3]
  )

  # Run model
  start_time = Sys.time()
  model.fit <-  R2jags::jags(
    model="/data/padgettn/diss_sim_results/model_code/dichotomous_items_cond_01_06.txt"
    , parameters.to.save = jags.params
    , data=mydata
    , n.chains = 4
    , n.burnin = 5000
    , n.iter = 10000
    , n.thin = 5
    , progress.bar = "none"
  )

  # extract summary information
  sample.summary <- model.fit$BUGSoutput$summary
  DIC <- model.fit$BUGSoutput$DIC
  pD <- model.fit$BUGSoutput$pD

  end_time = Sys.time()
  # edit Sys.time output for later
  elapsed_time <- difftime(end_time, start_time, units="mins")
  start_time <- as.character(start_time)
  end_time <- as.character(end_time)

  save.sys.info = matrix(c(c, i, start_time, end_time, elapsed_time, DIC, pD), nrow=1)
  save.samples <- as.data.frame(model.fit$BUGSoutput$sims.matrix)

  # check for output directory
  out.dir <- paste0("/data/padgettn/diss_sim_results/cond_",c,"/")
  if(dir.exists(out.dir) == F){
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  }
  # save posterior draws for one iteration (files are big)
  if(i == 1){
    suppressWarnings(write.table(save.samples, paste0(out.dir, "sim_posterior_draws_cond_",c,"_iter_",i,".txt"),
                                 sep="\t", col.names = T, row.names = T))
  }
  # Save posterior summary M, SD, quantiles
  suppressWarnings(write.table(sample.summary, paste0(out.dir, "sim_summary_cond_",c,"_iter_",i,".txt"),
                               sep="\t",col.names = T, row.names = T))
  # save model estimation information
  suppressWarnings(write.table(save.sys.info, paste0(out.dir, "sim_run_info_cond_",c,".txt"),
                               sep="\t", append=T, col.names = F))

  return(save.sys.info)

}


# thresholds
set.seed(1)
sim_tau <- matrix(
  runif(N_items, -0.5, 0.5),
  ncol=N_cat - 1, nrow=N_items, byrow=T
)

# Set up Condition 1
CON <- rep(list(list(c=COND_NUM,i=1, fx_sim_dat=simulate_data_misclass, N_persons=N_persons, N_items=N_items, N_cat = N_cat, tau=sim_tau)), N_rep)

for(i in 1:N_rep){
  CON[[i]]$i = i
}

# set number of parallel processes
registerDoParallel(36) # max number of ppn

# run model
system.time({
  RESULTS2 <- foreach(
    c=1:length(CON), .verbose = T, .combine = rbind) %dopar% ({
      ## Run the simulation for this row of the grid
      call_fx(CON[[c]])
    }) ## End foreach
}) ## End system.time
doParallel::stopImplicitCluster()

suppressWarnings(write.csv(RESULTS2, "/data/padgettn/diss_sim_results/cond_",COND_NUM,"/run_summary_cond_",COND_NUM,".csv", col.names = T, row.names = F))
q(save = "no")


