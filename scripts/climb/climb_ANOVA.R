############################################################
## Phenotype Analysis: Climbing Assay with Mixed Models
## Author: Leah Darwin
## Date: 2025-09-29
## Description: 
##   - Loads pre-processed climbing assay data
##   - Fits linear mixed-effects models (LMMs) using lmerTest
##   - Runs Type III ANOVAs
##   - Computes R² and effect sizes
##   - Generates publication-ready tables
############################################################

## ---------------------------
## Load required packages
## ---------------------------
packages <- c(
  "dplyr", "ggplot2", "lme4", "lmerTest", "broom", 
  "knitr", "kableExtra", "r2glmm", "performance"
)

# Install any missing packages
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------
## Global options
## ---------------------------
# Ensure Type III ANOVA with sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))

## ---------------------------
## Load data
## ---------------------------
climb <- read.csv("data/climb/climb_adj.csv", header = TRUE) %>%
  filter(!Mito %in% c("375", "Ore")) %>%   # Remove unwanted mito haplotypes
  mutate(
    across(c(Mito, Treatment, Nuc, Sex, Build, Vial), as.factor)
  )
  
## ---------------------------
## Stratify dataset by sex
## ---------------------------
climbF <- climb %>%
  filter(Sex == "F") %>%
  mutate(MitoNucBuild = paste(Mito, Nuc, Build, sep = ":"))

climbM <- climb %>%
  filter(Sex == "M") %>%
  mutate(MitoNucBuild = paste(Mito, Nuc, Build, sep = ":"))


## ---------------------------
## Mixed model: full dataset 
## ---------------------------
climb_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) + 
    Sex + Mito:Nuc:Sex + Sex:Treatment + Sex:Mito + Sex:Nuc + 
    (1 | Mito:Nuc:Treatment:Build:Vial),
  data = climb
)

# Type III ANOVA
climb_aov <- anova(climb_lm)

# Model fit R²
r2_nakagawa(climb_lm)

# Effect size R²
r2_climb <- r2beta(climb_lm)

# Helper: publication-ready kable for ANOVA tables
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

# Format ANOVA table with effect sizes
climb_aov_tab <- tidy(climb_aov) %>%
  select(-DenDF) %>%
  left_join(r2_climb %>% select(c(Effect, Rsq)), join_by(term == Effect)) %>%
  mutate(trait = "Climb") %>%
  select(trait, everything())

fmt_kable(climb_aov_tab)

## ---------------------------
## Female-only model
## ---------------------------
climbF_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) + 
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbF
)

climbF_aov <- anova(climbF_lm)
r2_nakagawa(climbF_lm)
r2_climbF <- r2beta(climbF_lm)

climbF_aov_tab <- tidy(climbF_aov) %>%
  select(-DenDF) %>%
  left_join(r2_climbF %>% select(c(Effect, Rsq)), join_by(term == Effect)) %>%
  mutate(trait = "Climbing") %>%
  select(trait, everything())

fmt_kable(climbF_aov_tab)

## ---------------------------
## Male-only model
## ---------------------------
climbM_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) + 
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbM
)

climbM_aov <- anova(climbM_lm)
r2_nakagawa(climbM_lm)
r2_climbM <- r2beta(climbM_lm)

climbM_aov_tab <- tidy(climbM_aov) %>%
  select(-DenDF) %>%
  left_join(r2_climbM %>% select(c(Effect, Rsq)), join_by(term == Effect)) %>%
  mutate(trait = "Climbing") %>%
  select(trait, everything())

fmt_kable(climbM_aov_tab)

