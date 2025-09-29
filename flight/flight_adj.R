library(dplyr)
library(ggplot2)

flight = read.csv("flight.csv", header=TRUE)

phenos_correct = flight %>% 
  na.omit() %>%
  group_by(Set,Sex,Treatment) %>%
  filter(Mito %in% c("yak","sm21") & Build=="A" | Build=="parental") %>%
  summarize(Set_adj = mean(Y))

##Get across set correction term
pos_correct = flight %>% 
  na.omit() %>%
  group_by(Treatment,Sex) %>%
  summarise(Pos_adj = mean(Y))

##Apply correction terms and filter out parentals and set0
flight_adj =  flight %>% 
  inner_join(phenos_correct, join_by(Set,Sex,Treatment)) %>%
  inner_join(pos_correct,join_by(Treatment,Sex)) %>%
  mutate(Y_adj = (Y-Set_adj)+Pos_adj)

plot_flight = flight_adj %>%
  group_by(Mito,Nuc,Set,Sex,Treatment,Build) %>%
  summarise(Y = mean(Y_adj))

hist(plot_flight$Y)
  
write.csv(flight_adj,"flight_adj.csv",row.names = FALSE,quote = FALSE)
