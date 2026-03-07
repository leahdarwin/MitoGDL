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
## Helper: publication-ready kable for ANOVA tables
## ---------------------------------------------------------
fmt_kable <- function(tab) {
  tab <- tab %>%
    mutate(p.value = format.pval(p.value, digits = 3, eps = 0.001))
  kable(
    tab,
    format    = "latex",
    booktabs  = TRUE,
    linesep   = "",
    escape    = FALSE,
    digits    = rep(4, 8),
    col.names = c("Trait", "Term", "SS", "MS", "$df$", "$F$", "$p$", "$R^2$")
  ) %>%
    kable_styling(
      latex_options = c("hold_position"),
      full_width    = FALSE,
      font_size     = 10
    ) %>%
    column_spec(2, width = "10em")
}

## ---------------------------------------------------------
## Format ANOVA results
## - Join with effect sizes
## - Clean column structure
## ---------------------------------------------------------
dev_aov <- tidy(dev_aov) %>%
  select(-DenDF) %>%
  left_join(r2_dev %>% select(c(Effect, Rsq)), join_by(term == Effect)) %>%
  mutate(trait = "Development") %>%
  select(trait, everything())

## ---------------------------------------------------------
## Export results table
## ---------------------------------------------------------
fmt_kable(dev_aov)
