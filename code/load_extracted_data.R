# easily grab data for all pages
mydata <- readr::read_csv(
  file = paste0(res_folder, "/simulation_summary_data_2022_02_12.csv"),
  col_types = "fiffffddddddddiddTTdfd")


# modify sigma_ts for N_cat = 2 conditions

for(i in 1:nrow(mydata)){
  if(is.na(mydata$parameter_group[i]) == F){
    if(mydata$parameter_group[i] == "sigma_st"){
      if(mydata$N_cat[i]!=2){
        mydata$true_value[i] = -mydata$true_value[i]
      }
    }
  }
}
