############################################################
# Script: weight_ANOVA.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Fit linear mixed models to adjusted fly weight 
#          data, evaluate fixed effects, compute effect sizes,
#          and format results for LaTeX-compatible tables.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("dplyr", "ggplot2", "lme4", "lmerTest", "broom",
              "knitr", "kableExtra", "r2glmm", "performance")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)


## ---------------------------------------------------------
## Set working options
## ---------------------------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))  # treatment contrasts for Type III ANOVA

## ---------------------------------------------------------
## Load and preprocess phenotype data
## ---------------------------------------------------------
weight <- read.csv("data/weight_adj.csv") %>%
  filter(!Mito %in% c("Ore", "375")) %>%  # remove unwanted mitochondrial lines
  mutate(Y_adj = Y_adj * 1000) # rescale phenotype
 

## Stratify by sex
weightF <- weight %>% filter(Sex == "F")
weightM <- weight %>% filter(Sex == "M")

## ---------------------------------------------------------
## Full model with sex interactions
## ---------------------------------------------------------
weight_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + 
    (1 | Mito:Nuc:Build) + 
    Sex + Mito:Nuc:Sex + Sex:Treatment + Sex:Mito + Sex:Nuc,
  data = weight
)

weight_aov <- anova(weight_lm)       # Type III ANOVA
weight_aov

r2_nakagawa(weight_lm)                # Marginal & conditional R²
r2_weight <- r2beta(weight_lm)       # Effect sizes (R² per term)

## Format p-values: fixed 4 decimal places for p >= 0.0001, switching to
## 2-significant-figure scientific notation below that (rather than the
## 0.0000 that fixed-decimal formatting would otherwise show).
format_p <- function(p, threshold = 1e-4) {
  ifelse(p < threshold,
         formatC(signif(p, 2), format = "e", digits = 1),
         formatC(p, format = "f", digits = 4))
}

## Helper: publication-ready kable for ANOVA tables
fmt_kable <- function(tab) {
  tab <- tab %>%
    mutate(p.value = format_p(p.value))
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

weight_aov <- tidy(weight_aov) %>%
  select(-DenDF) %>%
  left_join(r2_weight %>% select(Effect, Rsq), join_by(term == Effect)) %>%
  mutate(trait = "Weight") %>%
  select(trait, everything())

fmt_kable(weight_aov)

## ---------------------------------------------------------
## Female-only model
## ---------------------------------------------------------
weightF_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = weightF
)

weightF_aov <- anova(weightF_lm)
weightF_aov

r2_nakagawa(weightF_lm)
r2_weightF <- r2beta(weightF_lm)

weightF_aov <- tidy(weightF_aov) %>%
  select(-DenDF) %>%
  left_join(r2_weightF %>% select(Effect, Rsq), join_by(term == Effect)) %>%
  mutate(trait = "Weight") %>%
  select(trait, everything())

fmt_kable(weightF_aov)

## ---------------------------------------------------------
## Male-only model
## ---------------------------------------------------------
weightM_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment  + (1 | Mito:Nuc:Build),
  data = weightM
)

weightM_aov <- anova(weightM_lm)
weightM_aov

r2_nakagawa(weightM_lm)
r2_weightM <- r2beta(weightM_lm)

weightM_aov <- tidy(weightM_aov) %>%
  select(-DenDF) %>%
  left_join(r2_weightM %>% select(Effect, Rsq), join_by(term == Effect)) %>%
  mutate(trait = "Weight") %>%
  select(trait, everything())

fmt_kable(weightM_aov)
