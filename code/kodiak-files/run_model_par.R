#run_model_par

# this is a test
library(doParallel)
library(rjags)
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
    N = 500, J = 5, C = 5,
    etaCor = .23, etasd1 = 1, etasd2 = sqrt(0.1),
    lambda=0.9, nu=1.5, sigma.ei=0.25, rho1=0.1)
  sTau <- matrix(
    c(-2.5, -0.5, 1, 1.5,
      -2, -0.5, 1.5, 2.25,
      -2.25, -0.5, 1.5, 2.1,
      -1.75, -0.5, 2, 2.5,
      -1.75, -0.5, 1.75, 2.5),
    ncol=paravec[3]-1, nrow=paravec[2], byrow=T
  )
  sim.data <- simulate_data_misclass(paravec, tau=sTau)

  jags.params <- c(
    # ability measurement model
    "tau", "lambda", "theta", "psi",
    # reliability
    "reli.omega"
  )

  mydata <- list(
    y = sim.data$y,
    N = nrow(sim.data$y),
    nit=ncol(sim.data$y)
  )
  par.set <- parallel.seeds("base::BaseRNG", 2)
  # initialize
  model.fit <- jags.model(
    file="/data/padgettn/model1_ifa_jags.txt",
    #file=paste0(getwd(),"/code/kodiak-files/model1_ifa_jags.txt"),
    data=mydata,
    inits=par.set,
    n.chains = 2
  )

  a = Sys.time()
  model.fit$recompile()
  samples <- coda.samples(
    model.fit,
    variable.names = jags.params,
    n.iter = 100,
    thin=1,
    quiet=T
  )
  b = Sys.time()

  save.sys.info = c(c, i, Sys.info()[c("nodename","machine", "login")], b-a)
  save.samples <- as.matrix(samples,iters = T, chains = T)

  #suppressWarnings(write.csv(save.samples, paste0("data/sim_run_results_",c,"_",i,".csv"), append=T, col.names = FALSE, row.names = F))
  #if(i == 1){
  #  suppressWarnings(write.csv(colnames(save.samples), paste0("data/sim_results_colnames_",c,"_",i,".csv")))
  #}
  #suppressWarnings(write.csv(save.sys.info, "data/sim_run_info.csv", append=T, col.names = FALSE, row.names = F))

  suppressWarnings(write.csv(save.samples, paste0("/data/padgettn/sim_run_results_",c,"_",i,".csv"), append=T, col.names = FALSE, row.names = F))
  if(i==1){
    suppressWarnings(write.csv(colnames(save.samples), paste0("/data/padgettn/sim_results_colnames_",c,"_",i,".csv")))
  }
  suppressWarnings(write.csv(save.sys.info, "/data/padgettn/sim_run_info.csv", append=T, col.names = FALSE, row.names = F))

  return(save.sys.info)

}


CON <- list(
  list(c=1,i=1, fx_sim_dat=simulate_data_misclass),
  list(c=1,i=2, fx_sim_dat=simulate_data_misclass),
  list(c=1,i=3, fx_sim_dat=simulate_data_misclass),
  list(c=1,i=4, fx_sim_dat=simulate_data_misclass)
)


registerDoParallel(4) ## Use 4 cores
system.time({
  RESULTS2 <- foreach(
    c=1:length(CON), .verbose = T,.combine = rbind,
    .packages = c('rjags')) %dopar% ({
      ## Run the simulation for this row of the grid
      call_fx(CON[[c]])
    }) ## End foreach
}) ## End system.time
doParallel::stopImplicitCluster()

q(save = "no")


