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
options(contrasts = c("contr.sum", "contr.poly"))

#clade = read.csv("../mtDNAclade.csv")

##load phenotype 
#dist_df = read.csv("375Ore_distance.csv")

##load phenotype dataframes
weight = read.csv("weight/weight_adj.csv") %>%
  filter(!Mito%in%c("Ore","375")) %>%
  mutate(Y_adj = Y_adj*1000)

##Stratify on sex
weightF = weight %>% filter(Sex=="F")
weightM = weight %>% filter(Sex=="M") 

weight_lm = lmerTest::lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Mito:Nuc:Sex+Sex:Treatment+Sex:Mito+Sex:Nuc, data=weight)
weight_aov = anova(weight_lm)
weight_aov

r2_nakagawa(weight_lm)

r2_weight = r2beta(weight_lm)

weight_aov = tidy(weight_aov) %>%
  select(-DenDF) %>%
  left_join(r2_weight%>%select(c(Effect,Rsq)), join_by(term==Effect))

kable(weight_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

#weightC_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Clade:Nuc:Sex+Sex:Treatment+Sex:Clade+Sex:Nuc, data=weight)
#weightC_aov = anova(weightC_lm)
#weightC_aov
#weightC_aov = weightC_aov %>% select(c("NumDF", "Mean Sq", "F value", "Pr(>F)"))
#xtable(weightC_aov, digits=-2)

#weightDist_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Distance:Nuc:Sex+Sex:Treatment+Sex:Distance+Sex:Nuc, data=weight)
#weightDist_aov = anova(weightDist_lm)
#weightDist_aov

weightF_lm = lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build), data=weightF)
weightF_aov = anova(weightF_lm)
weightF_aov

r2_nakagawa(weightF_lm)

r2_weightF = r2beta(weightF_lm)

weightF_aov = tidy(weightF_aov) %>%
  select(-DenDF)%>%
  left_join(r2_weightF%>%select(c(Effect,Rsq)), join_by(term==Effect))

kable(weightF_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))



#weightFC_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build), data=weightF)
#weightFC_aov = anova(weightFC_lm)
#weightFC_aov

#weightFD_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build), data=weightF)
#weightFD_aov = anova(weightFD_lm)
#weightFD_aov

weightM_lm = lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build), data=weightM)
weightM_aov = anova(weightM_lm)
weightM_aov
r2_nakagawa(weightM_lm)

r2_weightM = r2beta(weightM_lm)

weightM_aov = tidy(weightM_aov) %>%
  select(-DenDF)%>%
  left_join(r2_weightM%>%select(c(Effect,Rsq)), join_by(term==Effect))
kable(weightM_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))


# weightMC_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build), data=weightM)
# weightMC_aov = anova(weightMC_lm)
# weightMC_aov
# 
# weightMD_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build), data=weightM)
# weightMD_aov = anova(weightMD_lm)
# weightMD_aov
# 
# weightF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   ggplot(aes(x = Distance, y = Y_adj, color=Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# weightM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   ggplot(aes(x = Distance, y = Y_adj, color=Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# weightFR = weightF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   filter(Treatment=="Rotenone")
# 
# weightF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   inner_join(weightFR, join_by(Mito,Nuc,Build,Distance)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y) %>%
#   ggplot(aes(x = Distance, y = CminusR, color=Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# weightMR = weightM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   filter(Treatment=="Rotenone")
# 
# weightM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   inner_join(weightFR, join_by(Mito,Nuc,Build,Distance)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y) %>%
#   ggplot(aes(x = Distance, y = CminusR, color=Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# treats = unique(weight$Treatment)
# nucs = unique(weight$Nuc)
# sexes = unique(weight$Sex)
# 
# for(sex in sexes){
#   for(treat in treats){
#     for(nuc in nucs){
#       temp_df = weight %>%
#         filter(Nuc==nuc & Treatment==treat & Sex==sex)
#       print(paste(sex, treat, nuc))
#       print(cor.test(temp_df$Distance, temp_df$Y_adj))
#     }
#   }
# }
# 
