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
weight <- read.csv("weight/data/weight_adj.csv") %>%
  filter(!Mito %in% c("Ore", "375")) %>%  # remove unwanted mitochondrial lines
  mutate(Y_adj = Y_adj * 1000)            # rescale phenotype

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

weight_aov <- tidy(weight_aov) %>%
  select(-DenDF) %>%
  left_join(r2_weight %>% select(Effect, Rsq), join_by(term == Effect))

kable(weight_aov, format = "latex", booktabs = TRUE,
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))

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
  left_join(r2_weightF %>% select(Effect, Rsq), join_by(term == Effect))

kable(weightF_aov, format = "latex", booktabs = TRUE,
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))

## ---------------------------------------------------------
## Male-only model
## ---------------------------------------------------------
weightM_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = weightM
)

weightM_aov <- anova(weightM_lm)

r2_nakagawa(weightM_lm)
r2_weightM <- r2beta(weightM_lm)

weightM_aov <- tidy(weightM_aov) %>%
  select(-DenDF) %>%
  left_join(r2_weightM %>% select(Effect, Rsq), join_by(term == Effect))

kable(weightM_aov, format = "latex", booktabs = TRUE,
      digits = c(0, 4, 4, 0, 4, 4, 4)) %>%
  kable_styling(latex_options = c("hold_position"))
