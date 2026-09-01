library(dplyr)
library(lme4)
library(lmerTest)

dev = read.csv("data/development_adj.csv")

ld_lm <- lmerTest::lmer(
  Larval_density ~ Mito * Nuc * Treatment  + (1 | Mito:Nuc:Build),
  data = dev
)
anova(ld_lm)
