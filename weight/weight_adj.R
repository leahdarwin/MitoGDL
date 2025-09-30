library(dplyr)
library(ggplot2)

weight = read.csv("weight.csv", header=TRUE)

##count the observations per group
counts = weight %>%
  group_by(Mito,Nuc,Treatment,Sex) %>%
  summarise(count=n(), Y = mean(weight_per_fly))


weight = weight %>%
  group_by(Sex) %>%
  mutate(z = scale(weight_per_fly)) %>%
  filter(abs(z) <= 3) %>%
  ungroup()

hist(weight$weight_per_fly)

##Get within set correction term
phenos_correct = weight %>% 
  na.omit() %>%
  group_by(Set,Sex,Treatment) %>%
  filter(Mito %in% c("yak","sm21") & Build=="A" | Build=="parental") %>%
  summarize(Set_adj = mean(weight_per_fly))

##Get across set correction term
pos_correct = weight %>% 
  na.omit() %>%
  group_by(Treatment,Sex) %>%
  summarise(Pos_adj = mean(weight_per_fly))

##Apply correction terms and filter out parentals and set0
weight_adj =  weight %>% 
  inner_join(phenos_correct, join_by(Set,Sex,Treatment)) %>%
  inner_join(pos_correct,join_by(Treatment,Sex)) %>%
  mutate(Y_adj = (weight_per_fly-Set_adj)+Pos_adj) 

write.csv(weight_adj,"weight_adj.csv",row.names = FALSE,quote = FALSE)

