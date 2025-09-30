############################################################
# Script: flight_ANOVA.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Perform mixed-model ANOVAs on adjusted flight 
#          phenotypes, stratified by sex, and extract effect 
#          sizes and R² values for reporting.
############################################################

# -------------------------------------------------------------------------
# Load required libraries (install if missing)
# -------------------------------------------------------------------------
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

# -------------------------------------------------------------------------
# Load and preprocess phenotype data
# -------------------------------------------------------------------------
# - Remove parental mito lines ("375", "Ore")
# - Average adjusted phenotype (Y_adj) across replicates
flight = read.csv("flight/data/flight_adj.csv") %>%
  filter(!Mito %in% c("375", "Ore")) %>%
  group_by(Mito, Nuc, Sex, Build, Treatment) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop")

# Stratify dataset by sex
flightF = flight %>% 
  filter(Sex == "F") %>%
  mutate(MitoNucBuild = paste(Mito, Nuc, Build, sep = ":"))

flightM = flight %>% 
  filter(Sex == "M") %>%
  mutate(MitoNucBuild = paste(Mito, Nuc, Build, sep = ":"))

# -------------------------------------------------------------------------
# Mixed-model analysis: all data
# -------------------------------------------------------------------------
# Model: Y_adj ~ Mito × Nuc × Treatment + Sex + interactions
# Random effects: (1 | Mito:Nuc:Build)
flight_lm = lmerTest::lmer(
  Y_adj ~ Mito * Nuc * Treatment +
    (1 | Mito:Nuc:Build) +
    Sex + Mito:Nuc:Sex + Sex:Treatment + Sex:Mito + Sex:Nuc,
  data = flight
)

# ANOVA table
flight_aov = anova(flight_lm)

# Marginal and conditional R²
r2_nakagawa(flight_lm)

# Effect sizes
r2_flight = r2beta(flight_lm)

# Merge ANOVA and effect size tables
flight_aov = tidy(flight_aov) %>%
  select(-DenDF) %>%
  left_join(r2_flight %>% select(c(Effect, Rsq)), join_by(term == Effect))

# Format for LaTeX output
kable(flight_aov, format = "latex", booktabs = TRUE, 
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))

# -------------------------------------------------------------------------
# Mixed-model analysis: females
# -------------------------------------------------------------------------
# Model with random effects: (1 | Mito:Nuc:Build)
flightF_lm = lmer(Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build), 
                  data = flightF)

# ANOVA table
flightF_aov = anova(flightF_lm)

# Marginal and conditional R²
r2_nakagawa(flightF_lm)

# Effect sizes
r2_flightF = r2beta(flightF_lm)

# Merge ANOVA and effect size tables
flightF_aov = tidy(flightF_aov) %>%
  select(-DenDF) %>%
  left_join(r2_flightF %>% select(c(Effect, Rsq)), join_by(term == Effect))

# Format for LaTeX output
kable(flightF_aov, format = "latex", booktabs = TRUE, 
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))

# -------------------------------------------------------------------------
# Mixed-model analysis: males
# -------------------------------------------------------------------------
# Model with random effects: (1 | Mito:Nuc:Build)
flightM_lm = lmer(Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build), 
                  data = flightM)

# ANOVA table
flightM_aov = anova(flightM_lm)

# Marginal and conditional R²
r2_nakagawa(flightM_lm)

# Effect sizes
r2_flightM = r2beta(flightM_lm)

# Merge ANOVA and effect size tables
flightM_aov = tidy(flightM_aov) %>%
  select(-DenDF) %>%
  left_join(r2_flightM %>% select(c(Effect, Rsq)), join_by(term == Effect))

# Format for LaTeX output
kable(flightM_aov, format = "latex", booktabs = TRUE, 
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))

