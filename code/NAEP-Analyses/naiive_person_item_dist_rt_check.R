# naiive_person_iem_dist_rt_check.R
# ===================================== #
# Math Motivation Study
# ===================================== #
# Created 
#   by: R. Noah Padgett
#   on: 2021-08-01
#
# Last Editted
#   by: R. Noah Padgett
#   on: 2021-08-06
# ===================================== #
# Purpose: 
# A crude approximation of person-item 
#   distance and response time 
#   relationship. Use CTT estimates.
#
# ===================================== #
# Load data
source("scripts/get_data.R")
# utility functions
source("scripts/get_utility_functions.R")

# ===================================== #
# Part 1. Use data subset (1% SRS)
N <- nrow(as.matrix(sdat_full[,mathIdent]))
person <- rowMeans(as.matrix(sdat_full[,mathIdent]),na.rm = T)
item <- colMeans(as.matrix(sdat_full[,mathIdent]),na.rm = T)
lrt <- pseudo_log(as.matrix(sdat_full[,colnames(sdat_full) %like% "RT_"]))
pidiff = person - item[1]

plot.dat <- data.frame(
  id = 1:N,
  person,
  item = paste0("item_",1),
  PIdiff = person - item[1],
  lrt = lrt[,1]
)
for(i in 2:5){
  p1 <-  data.frame(
    id = 1:N,
    person,
    item = paste0("item_",i),
    PIdiff = person - item[i],
    lrt = lrt[,i]
  )
  plot.dat <- full_join(plot.dat, p1)
}

cols = c("Parabolic Fit"="red", "GAM Fit"="blue")
p <- ggplot(plot.dat, aes(x=PIdiff, lrt, group=item))+
  geom_point()+
  geom_smooth(se=F,formula = y ~ poly(x,2), aes(color="Parabolic Fit"))+
  geom_smooth(se=F,aes(color='GAM Fit'))+
  facet_wrap(.~item, nrow=1)+
  scale_color_manual(values=cols, name="Fitted Regression Line")+
  labs(x="Person-Item Distance", y="log Response Time")+
  theme_classic()+
  theme(
    legend.position = c(0.8, 0.2)
  )
p
ggsave(plot=p, filename = "output/naiive_model_check_s.pdf",width = 7, heigh=5, units="in")

# 2D density plot 
p <- plot.dat %>% na.omit() %>%
  ggplot(aes(x=PIdiff, lrt, group=item))+
  geom_point(alpha=0.1, color="blue")+
  geom_density_2d(adjust=c(1.5,1.5), color="black")+
  facet_wrap(.~item, nrow=1)+
  labs(x="Person-Item Distance", y="log Response Time")+
  theme_classic()
p
ggsave(plot=p, filename = "output/naiive_model_check_dens_s.pdf",width = 7, heigh=5, units="in")

# ===================================== #
# Part 2. Use FULL sample of compete cases.
mydata_full_c <- mydataFull_RT_SurveyData_Combined %>%
  filter(
    accessionNumber %in% c("VH269048", "VH267478", "VH271749", "VH271337") |
      itemAccNum %in% c("VH268946")
  )
mydata_full_c_RT <- mydata_full_c %>%
  filter(itemAccNum %in% mathIdentV) %>%
  mutate(itemAccRecode = factor(itemAccNum, levels = mathIdentV, labels=mathIdent)) %>%
  select(studentIDnew,  itemAccRecode, RT, Revisit) %>%
  group_by(studentIDnew,  itemAccRecode, .drop=F) %>%
  summarise(RT = sum(RT)) %>%
  pivot_wider(
    id_cols = c("studentIDnew"),
    names_from=itemAccRecode,
    values_from = RT,
    names_prefix = "RT_" 
  )
mydata_full_c_IR = mydata_full_c %>%
  filter(itemAccNum %in% mathIdentV) %>%
  select(studentIDnew, all_of(mathIdent)) %>%
  unique()
# combine data into one object to make sure the IDs align
mydata_full_c <- full_join(mydata_full_c_RT, mydata_full_c_IR) %>%
  na.omit() # remove cases with missing values - complicates things


mydata_full_c <- mydata_full_c %>%
  select(all_of(mathIdent), all_of(grep("RT_", colnames(mydata_full_c),value = T))) %>%
  na.omit()
N <- nrow(as.matrix(mydata_full_c[,mathIdent]))
person <- rowMeans(as.matrix(mydata_full_c[,mathIdent]),na.rm = T)
item <- colMeans(as.matrix(mydata_full_c[,mathIdent]),na.rm = T)
lrt <- pseudo_log(as.matrix(mydata_full_c[,colnames(mydata_full_c) %like% "RT_"]))
pidiff = person - item[1]

plot.dat <- data.frame(
  id = 1:N,
  person,
  item = paste0("item_",1),
  PIdiff = person - item[1],
  lrt = lrt[,1]
)
for(i in 2:5){
  p1 <-  data.frame(
    id = 1:N,
    person,
    item = paste0("item_",i),
    PIdiff = person - item[i],
    lrt = lrt[,i]
  )
  plot.dat <- full_join(plot.dat, p1)
}

cols = c("Parabolic Fit"="red", "GAM Fit"="blue")
p <- ggplot(plot.dat, aes(x=PIdiff, lrt, group=item))+
  geom_smooth(se=F,formula = y ~ poly(x,2), aes(color="Parabolic Fit"))+
  geom_smooth(se=F,aes(color='GAM Fit'))+
  facet_wrap(.~item, nrow=1)+
  scale_color_manual(values=cols, name="Fitted Regression Line")+
  labs(x="Person-Item Distance", y="log Response Time")+
  theme_classic()+
  theme(
    legend.position = c(0.8, 0.2)
  )
p
ggsave(plot=p, filename = "output/naiive_model_check_f.pdf",width = 7, heigh=5, units="in")

# 2D density plot
p <- ggplot(plot.dat, aes(x=PIdiff, lrt, group=item))+
  geom_point(alpha=0.1, color="blue")+
  geom_density_2d(adjust=c(1.5,1.5), color="black")+
  facet_wrap(.~item, nrow=1)+
  labs(x="Person-Item Distance", y="log Response Time")+
  theme_classic()
p
ggsave(plot=p, filename = "output/naiive_model_check_dens_f.pdf",width = 7, heigh=5, units="in")
