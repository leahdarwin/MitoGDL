library(dplyr)
library(ggplot2)

##load development csv
dev = read.csv("development.csv")
dev = na.omit(dev_dat)
hist(dev$TSE)

##Get within set correction term
phenos_correct = dev %>% 
  na.omit() %>%
  group_by(Set,Treatment) %>%
  filter(Mito %in% c("yak","sm21") & Build=="A" | Build=="parental") %>%
  summarize(Set_adj = mean(TSE))

##Get across set correction term
pos_correct = dev %>% 
  na.omit() %>%
  group_by(Treatment) %>%
  summarise(Pos_adj = mean(TSE))

total_treat = dev %>%
  na.omit() %>%
  group_by(Treatment) %>%
  summarise(Total_treat = mean(Total))

##Apply correction terms and filter out parentals and set0
dev_adj =  dev %>% 
  inner_join(phenos_correct, join_by(Set,Treatment)) %>%
  inner_join(pos_correct,join_by(Treatment)) %>%
  inner_join(total_treat, join_by(Treatment)) %>%
  mutate(Y_adj = (TSE-Set_adj)+Pos_adj) %>%
  mutate(Larval_density = Total/Total_treat)

hist(dev_adj$Y_adj)

write.csv(dev_adj, "development_adj.csv", row.names=FALSE, quote=FALSE)
