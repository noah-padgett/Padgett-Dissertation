# Dissertation simulation summary
library(tidyverse)
library(data.table)
options(scipen = 999, digits=3)

# results folder
res_folder <- "C:/Users/noahp/Desktop/dissertation_sim_results"

# CONDITIONS
CONDITIONS <- list(
  N = rep(c(500, 2500), 6),
  N_cat = c(rep(2,6), rep(5,6)), # , rep(3,6) , rep(7,6)
  N_items = rep(c(5,5,10,10,20,20),2)
)

# extract info data.frame
Ncon = 12
Niter = 5
Npara =  20*4+20*4+4 # just need max
nSkip = c(0, 36, 66, 66, 126, 126, 45, 45, 85, 85, 165, 165)
cNames <- c(
  "condition", "iter", "N", "N_items", "N_cat", "parameter", "post_mean", "post_sd", "post_q025", "post_q25", "post_q50", "post_q75", "post_q975", "Rhat", "n.eff", "DIC", "pD", "start_time", "end_time", "elapsed_time"
)
res_out <- as.data.frame(matrix(ncol=length(cNames), nrow=Ncon*Niter*Npara)) # will be too big but easy to reduce later filled it
colnames(res_out) <- cNames
cond_number <- 1
iter <- 1
res_dat_row_iterator <- 1

for(cond_number in 1:Ncon){
  # extract condition info
  N <- CONDITIONS$N[cond_number]
  N_cat <- CONDITIONS$N_cat[cond_number]
  N_items <- CONDITIONS$N_items[cond_number]

  # get run infor to fill in
  run_info <- read.delim(paste0(res_folder, "/cond_",cond_number,"/sim_run_info_cond_", cond_number,".txt"), header=FALSE)
  colnames(run_info) <- c("skip", "condition", "iteration", "start_time", "end_time", "elapsed_time", "DIC", "pD")

  for(iter in 1:Niter){

    # subset to iteration
    run_info_i <- filter(run_info, iteration == iter)

    # get posterior summary data
    if(iter <= 5){
      post_summary <- read.delim(paste0(res_folder, "/cond_",cond_number,"/sim_summary_cond_", cond_number,"_iter_",iter,".txt"), row.names=1,skip = nSkip[cond_number])
    } else {
      post_summary <- read.delim(paste0(res_folder, "/cond_",cond_number,"/sim_summary_cond_", cond_number,"_iter_",iter,".txt"), row.names=1)
    }

    post_sum_nrow <- nrow(post_summary)
    post_sum_rownames <- rownames(post_summary)

    # fill in results dataframe
    iter_ps <- 1
    for(iter_ps in 1:post_sum_nrow){
      res_out[res_dat_row_iterator, 1] <- cond_number
      res_out[res_dat_row_iterator, 2] <- iter
      res_out[res_dat_row_iterator, 3] <- N
      res_out[res_dat_row_iterator, 4] <- N_items
      res_out[res_dat_row_iterator, 5] <- N_cat
      res_out[res_dat_row_iterator, 6] <- post_sum_rownames[iter_ps]
      res_out[res_dat_row_iterator, 7:15] <- post_summary[iter_ps, ]
      res_out[res_dat_row_iterator, 16:20] <- run_info_i[c("DIC", "pD", "start_time", "end_time", "elapsed_time")]

      # update row number
      res_dat_row_iterator = res_dat_row_iterator + 1
    }
  }
  cat("\n Condition:\t", cond_number)
}

res_out <- na.omit(res_out)


res_out %>%
  group_by(condition, parameter) %>%
  summarise(
    avg_post_M = mean(post_mean),
    avg_post_SD = mean(post_sd),
    avg_Rhat = mean(Rhat),
    avg_neff = mean(n.eff)
  ) %>%
  filter(parameter %like% "lambda", !(parameter %like% "lambda.std")) %>%
  View()


# thresholds
set.seed(1)
cond_number <- 1
iter <- 1
N <- CONDITIONS$N[cond_number]
N_cat <- CONDITIONS$N_cat[cond_number]
N_items <- CONDITIONS$N_items[cond_number]
sim_tau <- matrix(
  runif(N_items, -0.5, 0.5),
  ncol=N_cat - 1, nrow=N_items, byrow=T
)

