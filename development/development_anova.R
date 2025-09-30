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

# clade = read.csv("../mtDNAclade.csv")
# dist_df = read.csv("375Ore_distance.csv")
# distYak_df = read.csv("yak_distance.csv") %>% rename(DistY = Distance)

##load phenotype dataframes
dev = read.csv("development/development_adj.csv") %>%
  filter(!Mito%in%c("375","Ore"))

dev_lm = lmerTest::lmer(Y_adj ~ Mito*Nuc*Treatment+(1|Mito:Nuc:Build) +Larval_density, data=dev)
dev_aov = anova(dev_lm)
dev_aov

r2_nakagawa(dev_lm)

##get effect sizes 
r2_dev = r2beta(dev_lm)

dev_aov = tidy(dev_aov) %>%
  select(-DenDF)%>%
  left_join(r2_dev%>%select(c(Effect,Rsq)), join_by(term==Effect))

kable(dev_aov, format = "latex", booktabs = TRUE, digits=c(0,4,4,0,4,4,4)) %>%
  kable_styling(latex_options = c("hold_position"))

# devM_lm = lmerTest::lmer(Y_adj ~ Clade*Nuc*Treatment+(1|Mito:Nuc:Build) + (1|Mito:Nuc:Build:Treatment:Vial)+Larval_density, data=dev)
# devM_aov = anova(devM_lm)
# devM_aov
# 
# devDist_lm = lmerTest::lmer(Y_adj ~ Distance*Nuc*Treatment+(1|Mito:Nuc:Build) + (1|Mito:Nuc:Build:Treatment:Vial)+Larval_density, data=dev)
# devDist_aov = anova(devDist_lm)
# devDist_aov
# 
# devDistY_lm = lmerTest::lmer(Y_adj ~ DistY*Nuc*Treatment+(1|Mito:Nuc:Build) + (1|Mito:Nuc:Build:Treatment:Vial)+Larval_density, data=dev)
# devDistY_aov = anova(devDistY_lm)
# devDistY_aov
# 
# dev_R = dev %>%
#   group_by(Mito,Nuc,Treatment,Distance,Build) %>%
#   summarise(Y_adj = mean(Y_adj) )%>%
#   filter(Treatment=="Rotenone")
# dev_CminusR = dev %>%
#   group_by(Mito,Nuc,Treatment,Distance,Build) %>%
#   summarise(Y_adj = mean(Y_adj) )%>%
#   filter(Treatment=="Control") %>%
#   inner_join(dev_R, join_by(Mito,Nuc,Build)) %>%
#   mutate(CminusR = Y_adj.x - Y_adj.y)
# 
# dev %>%
#   group_by(Mito,Nuc,Build,Treatment,Distance, DistY) %>%
#   summarise(Y_adj = mean(Y_adj) )%>%
#  # ggplot(aes(x = Distance, y = Y_adj, colour = Treatment)) +
#   ggplot(aes(x = DistY, y = Y_adj, colour = Treatment)) +
#   facet_wrap(~Nuc) +
#   geom_point()+
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw() +
#   ylab("Development Time (days)")
# 
# dev_CminusR %>%
#   ggplot(aes(x=Distance.x, y = CminusR, color = Nuc)) +
#   geom_point() +
#   stat_poly_line() +
#   stat_poly_eq(use_label(c("eq", "R2"))) +
#   theme_linedraw() +
#   ylab("Treatment Effect (Control-Rotenone)")
# 
# treats = unique(dev$Treatment)
# nucs = unique(dev$Nuc)
# 
# for(treat in treats){
#   for(nuc in nucs){
#       temp_df = dev %>%
#         filter(Nuc==nuc & Treatment==treat)
#       print(paste(treat, nuc))
#       #print(cor.test(temp_df$Distance, temp_df$Y_adj))
#       print(cor.test(temp_df$DistY, temp_df$Y_adj))
#     }
# }
# 
