library(dplyr)
library(ggplot2)
library(ggfortify)
library(ggsci)
library(grid)
library(gridExtra)
library(pheatmap)
library(vegan) ##for procrustes 


setwd("/Users/leahdarwin/Documents/drand/gdlPhenotypeExperiment/dist")

weight = read.csv("weight/weight.csv")
climb = read.csv("climb/climb.csv")
dev = read.csv("development/development.csv")
flight = read.csv("flight/flight.csv")

weight_F = weight %>%
  na.omit() %>%
  filter(Sex == "F") %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(weight_F = mean(weight_per_fly))

weight_M = weight %>%
  na.omit() %>%
  filter(Sex == "M") %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(weight_M = mean(weight_per_fly))

climb_F = climb %>%
  na.omit() %>%
  filter(Sex == "F") %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(climb_F = mean(Slope))

climb_M = climb %>%
  na.omit() %>%
  filter(Sex == "M") %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(climb_M = mean(Slope))

flight_F = flight %>%
  na.omit() %>%
  filter(Sex == "F") %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(flight_F = mean(Y))

flight_M = flight %>%
  na.omit() %>%
  filter(Sex == "M") %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(flight_M = mean(Y))

dev_MF = dev %>%
  na.omit() %>%
  group_by(Mito,Nuc,Treatment) %>%
  summarise(dev = mean(TSE))

merged = full_join(dev_MF, full_join(weight_M, full_join(climb_F, full_join(climb_M, full_join(flight_F, full_join(flight_M, weight_F)))))) %>% na.omit()

merged = merged %>%
  mutate(Mito_orig = case_when(grepl("Bei", Mito) ~ "Bei",
                               grepl("Z", Mito) ~ "Zim",
                               .default = Mito))

merged_F_control = merged %>%
  filter(Treatment == "Control") %>%
  select(-c(weight_M, climb_M, flight_M)) %>%
  group_by(Mito) %>%
  summarise_if(is.numeric, mean)

merged_F_rote = merged %>%
  filter(Treatment == "Rotenone") %>%
  select(-c(weight_M, climb_M, flight_M)) %>%
  group_by(Mito) %>%
  summarise_if(is.numeric, mean)

merged_F_Ore = merged %>%
  filter(Nuc=="Ore") %>%
  select(-c(weight_M, climb_M, flight_M))

merged_F_375 = merged %>%
  filter(Nuc=="375") %>%
  select(-c(weight_M, climb_M, flight_M))

merged_F_control_nuc = merged %>%
  filter(Treatment == "Control") %>%
  select(-c(weight_M, climb_M, flight_M)) %>%
  group_by(Nuc) %>%
  summarise_if(is.numeric, mean)

merged_F_rote_nuc = merged %>%
  filter(Treatment == "Rotenone") %>%
  select(-c(weight_M, climb_M, flight_M)) %>%
  group_by(Nuc) %>%
  summarise_if(is.numeric, mean)

merged_F_rote_mitonuc = merged %>%
  filter(Treatment == "Rotenone") %>%
  select(-c(weight_M, climb_M, flight_M))

merged_F_control_mitonuc = merged %>%
  filter(Treatment == "Control") %>%
  select(-c(weight_M, climb_M, flight_M))

pca = prcomp(merged[,4:10], scale. = TRUE)

pca_F_control = prcomp(merged_F_control[,2:5], scale. = TRUE)
pca_F_rote = prcomp(merged_F_rote[,2:5], scale. = TRUE)
pca_F_control_mitonuc = prcomp(merged_F_control_mitonuc[,4:7], scale. = TRUE)
pca_F_rote_mitonuc = prcomp(merged_F_rote_mitonuc[,4:7], scale. = TRUE)
pca_F_Ore = prcomp(merged_F_Ore[,4:7], scale. = TRUE)
pca_F_375 = prcomp(merged_F_375[,4:7], scale. = TRUE)

autoplot(pca, data=merged, color = 'Nuc', x = 1, y = 2)+
  scale_color_manual(values = c("Ore" = "#67a9cf", "375"="#ef8a62")) +
  theme_light() +
  labs(color="Nuclear Genotype")+ 
  theme(legend.position="top")

autoplot(pca_F_control, color='Mito',data=merged_F_control, x=1, y=2)
autoplot(pca_F_rote, color='Mito',data=merged_F_rote, x=1, y=2)
  
pro_F_treat = procrustes(pca_F_control, pca_F_rote, symmetric = TRUE)
plot(pro_F_treat, kind=1)
plot(pro_F_treat, kind=2)
protest(pca_F_control, pca_F_rote)

pro_F_nuc = procrustes(pca_F_Ore, pca_F_375, symmetric = TRUE)
plot(pro_F_nuc, kind=1)
plot(pro_F_nuc, kind=2)
protest(pca_F_Ore, pca_F_375)

pro_F_treat_mitonuc = procrustes(pca_F_control_mitonuc, pca_F_rote_mitonuc, symmetric = TRUE)
plot(pro_F_treat_mitonuc, kind=1)
plot(pro_F_treat_mitonuc, kind=2)
protest(pca_F_control_mitonuc, pca_F_rote_mitonuc)
