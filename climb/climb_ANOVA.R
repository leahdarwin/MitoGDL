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

##set global effects coding for typeIII anovas
options(contrasts = c("contr.sum", "contr.poly"))

#clade = read.csv("../mtDNAclade.csv")

#dist_df = read.csv("375Ore_distance.csv")

##load phenotype dataframes
climb = read.csv("climb/climb_adj.csv") %>%
  filter(!Mito %in% c("375","Ore")) %>%
  mutate(Mito = as.factor(Mito)) %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(Nuc = as.factor(Nuc)) %>%
  mutate(Sex = as.factor(Sex)) %>%
  mutate(Build = as.factor(Build)) %>%
  mutate(Vial = as.factor(Vial))
  

##Stratify on sex
climbF = climb %>% filter(Sex=="F") %>% mutate(MitoNucBuild = paste(Mito,Nuc,Build,sep=":"))
climbM = climb %>% filter(Sex=="M") %>% mutate(MitoNucBuild = paste(Mito,Nuc,Build,sep=":"))


unique_mito = unique(climb$Mito)
subset_mito = sample(unique_mito, size = 10, replace = FALSE)
climb_subset = climb[climb$Mito %in% subset_mito, ]

climb_lm = lmerTest::lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Mito:Nuc:Sex+Sex:Treatment+Sex:Mito+Sex:Nuc+(1|Mito:Nuc:Treatment:Build:Vial), data=climb_subset)
climb_aov = anova(climb_lm)
climb_aov

##get marginal and conditional r2 
r2_nakagawa(climb_lm)

##get effect sizes 
r2_climb = r2beta(climb_lm)

##make table
climb_aov_tab = tidy(climb_aov) %>%
  select(-DenDF) %>%
  left_join(r2_climb%>%select(c(Effect,Rsq)), join_by(term==Effect))
kable(climb_aov_tab, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))


climbF_lm = lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build)+(1|Mito:Nuc:Build:Treatment:Vial), data=climbF)
climbF_aov = anova(climbF_lm)
climbF_aov

r2_nakagawa(climbF_lm)
r2_climbF = r2beta(climbF_lm)

climbF_aov = tidy(climbF_aov) %>%
  select(-DenDF) %>%
  left_join(r2_climbF%>%select(c(Effect,Rsq)), join_by(term==Effect))
kable(climbF_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

climbM_lm = lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build), data=climbM)
climbM_aov = anova(climbM_lm)
climbM_aov

r2_nakagawa(climbM_lm)
r2_climbM = r2beta(climbM_lm)

climbM_aov = tidy(climbM_aov) %>%
  select(-DenDF) %>%
  left_join(r2_climbM%>%select(c(Effect,Rsq)), join_by(term==Effect))
kable(climbM_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

# climbC_lm = lmerTest::lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Clade:Nuc:Sex+Sex:Treatment+Sex:Clade+Sex:Nuc + (1|Mito:Nuc:Build:Sex:Treatment:Vial), data=climb)
# climbC_aov = anova(climbC_lm)
# climbC_aov
# climbC_aov = climbC_aov %>% select(c("NumDF", "Mean Sq", "F value", "Pr(>F)"))
# xtable(climbC_aov, digits=-2)
# 
# climbCF_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build)+ (1|Mito:Nuc:Build:Sex:Treatment:Vial), data=climbF)
# climbCF_aov = anova(climbCF_lm)
# climbCF_aov
# climbCF_aov = climbCF_aov %>% select(c("NumDF", "Mean Sq", "F value", "Pr(>F)"))
# xtable(climbCF_aov, digits=-2)
# 
# climbCM_lm = lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build)+ (1|Mito:Nuc:Build:Sex:Treatment:Vial), data=climbM)
# climbCM_aov = anova(climbCM_lm)
# climbCM_aov
# climbCM_aov = climbCM_aov %>% select(c("NumDF", "Mean Sq", "F value", "Pr(>F)"))
# xtable(climbCM_aov, digits=-2)
# 
# climbDist_lm = lmerTest::lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build)+Sex+Distance:Nuc:Sex+Sex:Treatment+Sex:Distance+Sex:Nuc + (1|Mito:Nuc:Build:Sex:Treatment:Vial), data=climb)
# climbDist_aov = anova(climbDist_lm)
# climbDist_aov
# climbDist_aov = climbDist %>% select(c("NumDF", "Mean Sq", "F value", "Pr(>F)"))
# xtable(climbC_aov, digits=-2)
# 
# climbDistF_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build)+ (1|Mito:Nuc:Build:Sex:Treatment:Vial), data=climbF)
# climbDistF_aov = anova(climbDistF_lm)
# climbDistF_aov
# 
# climbDistM_lm = lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build)+ (1|Mito:Nuc:Build:Sex:Treatment:Vial), data=climbM)
# climbDistM_aov = anova(climbDistM_lm)
# climbDistM_aov
# 
# climbF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   ggplot(aes(x = Distance, y = Y_adj, color=Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# climbFR = climbF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   filter(Treatment=="Rotenone")
# 
# climbF %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   inner_join(climbFR, join_by(Mito,Nuc,Build,Distance)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y) %>%
#   ggplot(aes(x = Distance, y = CminusR, color=Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# 
# climbM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   ggplot(aes(x = Distance, y = Y_adj, color=Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# 
# climbMR = climbM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   filter(Treatment=="Rotenone")
# 
# climbM %>%
#   group_by(Mito,Nuc,Treatment,Build,Distance) %>%
#   summarise(Y_adj = mean(Y_adj))%>%
#   inner_join(climbMR, join_by(Mito,Nuc,Build,Distance)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y) %>%
#   ggplot(aes(x = Distance, y = CminusR, color=Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw()
# 
# treats = unique(climb$Treatment)
# nucs = unique(climb$Nuc)
# sexes = unique(climb$Sex)
# 
# for(sex in sexes){
#   for(treat in treats){
#     for(nuc in nucs){
#       temp_df = climb %>%
#         filter(Nuc==nuc & Treatment==treat & Sex==sex)
#       print(paste(sex, treat, nuc))
#       print(cor.test(temp_df$Distance, temp_df$Y_adj))
#     }
#   }
# }
