############################################################
# Script: development_anova.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Fit linear mixed models to adjusted pupation 
#          (development) data, evaluate fixed effects, 
#          compute effect sizes, and format results for 
#          LaTeX-compatible tables.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c(
  "dplyr", "ggplot2", "lme4", "lmerTest", 
  "broom", "knitr", "kableExtra", "r2glmm", "performance"
)
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------------------------------------
## Load adjusted development dataset
## - Exclude unwanted mitochondrial backgrounds (375, Ore)
## ---------------------------------------------------------
dev <- read.csv("development/data/development_adj.csv") %>%
  filter(!Mito %in% c("375", "Ore"))

## ---------------------------------------------------------
## Fit linear mixed-effects model
## - Fixed effects: Mito × Nuc × Treatment, Larval_density
## - Random effect: Mito:Nuc:Build (genotype nested in build)
## ---------------------------------------------------------
dev_lm <- lmerTest::lmer(
  Y_adj ~ Mito * Nuc * Treatment + Larval_density + (1 | Mito:Nuc:Build),
  data = dev
)

## ---------------------------------------------------------
## Extract ANOVA table
## ---------------------------------------------------------
dev_aov <- anova(dev_lm)

## ---------------------------------------------------------
## Compute marginal and conditional R²
## - Nakagawa’s method
## ---------------------------------------------------------
r2_nakagawa(dev_lm)

## ---------------------------------------------------------
## Compute effect sizes (semi-partial R² for fixed effects)
## ---------------------------------------------------------
r2_dev <- r2beta(dev_lm)

## ---------------------------------------------------------
## Format ANOVA results
## - Join with effect sizes
## - Clean column structure
## ---------------------------------------------------------
dev_aov <- tidy(dev_aov) %>%
  select(-DenDF) %>%
  left_join(r2_dev %>% select(c(Effect, Rsq)), join_by(term == Effect))

## ---------------------------------------------------------
## Export results table
## - LaTeX-compatible table with kable
## ---------------------------------------------------------
kable(dev_aov, format = "latex", booktabs = TRUE,
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))
