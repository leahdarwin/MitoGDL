library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom)
library(knitr)
library(kableExtra)
library(r2glmm)
library(performance)

setwd("/Users/leahdarwin/Documents/drand/gdlPhenotypeExperiment/dist")

#clade = read.csv("../mtDNAclade.csv")

##load phenotype 
#dist_df = read.csv("375Ore_distance.csv")

##load phenotype dataframes
flight = read.csv("flight/flight_adj.csv") %>%
  filter(!Mito %in% c("375","Ore")) %>%
  group_by(Mito,Nuc,Sex,Build,Treatment) %>%
  summarise(Y_adj = mean(Y_adj))

##Stratify on sex
flightF = flight %>% filter(Sex=="F") %>% mutate(MitoNucBuild = paste(Mito,Nuc,Build,sep=":"))
flightM = flight %>% filter(Sex=="M") %>% mutate(MitoNucBuild = paste(Mito,Nuc,Build,sep=":"))

flight_lm = lmerTest::lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Mito:Nuc:Sex+Sex:Treatment+Sex:Mito+Sex:Nuc, data=flight)
flight_aov = anova(flight_lm)
flight_aov

##get marginal and conditional r2 
r2_nakagawa(flight_lm)

##get effect sizes 
r2_flight = r2beta(flight_lm)

flight_aov = tidy(flight_aov) %>%
  select(-DenDF)%>%
  left_join(r2_flight%>%select(c(Effect,Rsq)), join_by(term==Effect))

kable(flight_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

# flightC_lm = lmerTest::lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Clade:Nuc:Sex+Sex:Treatment+Sex:Clade+Sex:Nuc, data=flight)
# flightC_aov = anova(flightC_lm)
# flightC_aov
# 
# flightD_lm = lmerTest::lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Distance:Nuc:Sex+Sex:Treatment+Sex:Distance+Sex:Nuc, data=flight)
# flightD_aov = anova(flightD_lm)
# flightD_aov

flightF_lm = lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build), data=flightF)
flightF_aov = anova(flightF_lm)
flightF_aov

r2_nakagawa(flightF_lm)

##get effect sizes 
r2_flightF = r2beta(flightF_lm)

flightF_aov = tidy(flightF_aov) %>%
  select(-DenDF)%>%
  left_join(r2_flightF%>%select(c(Effect,Rsq)), join_by(term==Effect))
kable(flightF_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

# flightFC_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build), data=flightF)
# flightFC_aov = anova(flightFC_lm)
# flightFC_aov
# 
# flightFD_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build), data=flightF)
# flightFD_aov = anova(flightFD_lm)
# flightFD_aov

flightM_lm = lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build), data=flightM)
flightM_aov = anova(flightM_lm)
flightM_aov

##get marginal and conditional r2 
r2_nakagawa(flightM_lm)

##get effect sizes 
r2_flightM = r2beta(flightM_lm)

flightM_aov = tidy(flightM_aov) %>%
  select(-DenDF)%>%
  left_join(r2_flightM%>%select(c(Effect,Rsq)), join_by(term==Effect))
kable(flightM_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

 

# flightMC_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build), data=flightM)
# flightMC_aov = anova(flightMC_lm)
# flightMC_aov
# 
# flightMD_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build), data=flightM)
# flightMD_aov = anova(flightMD_lm)
# flightMD_aov
# 
# flightF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   ggplot(aes(x = Distance, y = Y_adj, color=Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# flightM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   ggplot(aes(x = Distance, y = Y_adj, color=Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# flightFR = flightF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   filter(Treatment=="Rotenone")
# 
# flightF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   inner_join(flightFR, join_by(Mito,Nuc,Build,Distance)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y) %>%
#   ggplot(aes(x = Distance, y = CminusR, color=Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# flightMR = flightM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   filter(Treatment=="Rotenone")
# 
# flightM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   inner_join(flightMR, join_by(Mito,Nuc,Build,Distance)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y) %>%
#   ggplot(aes(x = Distance, y = CminusR, color=Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# treats = unique(flight$Treatment)
# nucs = unique(flight$Nuc)
# sexes = unique(flight$Sex)
# 
# for(sex in sexes){
#   for(treat in treats){
#     for(nuc in nucs){
#       temp_df = flight %>%
#         filter(Nuc==nuc & Treatment==treat & Sex==sex)
#       print(paste(sex, treat, nuc))
#       print(cor.test(temp_df$Distance, temp_df$Y_adj))
#     }
#   }
# }
