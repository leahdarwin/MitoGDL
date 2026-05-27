############################################################
# Script: surv_ANOVA.R
# Author: Leah Darwin
# Date: 2026-05-05
# Purpose: Perform ANOVAs on adult survival counts and produce
#          separate LaTeX ANOVA and estimated marginal means 
#          tables via kable booktabs.
############################################################

# -------------------------------------------------------------------------
# Load required libraries (install if missing)
# -------------------------------------------------------------------------
packages <- c("dplyr", "broom", "emmeans", "lme4", "lmerTest", "knitr", "kableExtra", "r2glmm", "performance")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p, dependencies = TRUE)
}
lapply(packages, library, character.only = TRUE)

# Ensure Type III ANOVA with sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------
egg_count <- read.csv("data/survival.csv") %>% mutate(Total = F_total + M_total)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

# Tidy ANOVA table from lmerTest model with R² effect sizes
tidy_aov <- function(model) {
  r2 <- r2beta(model) %>% select(Effect, Rsq)
  tidy(anova(model)) %>%
    select(term, sumsq, meansq, NumDF, statistic, p.value) %>%
    left_join(r2, by = join_by(term == Effect)) %>%
    mutate(
      across(c(sumsq, meansq, statistic, Rsq), ~round(., 4)),
      p.value = format.pval(p.value, digits = 3, eps = 0.001)
    )
}

# Tidy EMMs using the supplied spec; keeps all columns emmeans returns
tidy_emm <- function(model, spec) {
  as.data.frame(emmeans(model, spec)) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
}

# kable for ANOVA results, rows grouped by trait via pack_rows
anova_table <- function(models, labels) {
  aov_list <- lapply(models, tidy_aov)
  aov_n    <- sapply(aov_list, nrow)
  combined <- bind_rows(aov_list)

  tbl <- kable(
    combined,
    format    = "latex",
    booktabs  = TRUE,
    linesep   = "",
    escape    = FALSE,
    col.names = c("Term", "SS", "MS", "$df$", "$F$", "$p$", "$R^2$")
  ) %>%
    kable_styling(
      latex_options = c("hold_position"),
      full_width    = FALSE,
      font_size     = 10
    ) %>%
    column_spec(1, width = "10em")

  cur <- 1
  for (i in seq_along(labels)) {
    tbl <- pack_rows(tbl, labels[i], cur, cur + aov_n[i] - 1)
    cur <- cur + aov_n[i]
  }
  tbl
}

# kable for EMMs, rows grouped by trait via pack_rows
# spec: emmeans formula, e.g. ~ Treatment or ~ Treatment + Egg_Count
# col.names: must match the columns returned by the chosen spec
emm_table <- function(models, labels, spec, col.names) {
  emm_list <- lapply(models, tidy_emm, spec = spec)
  emm_n    <- sapply(emm_list, nrow)
  combined <- bind_rows(emm_list)

  tbl <- kable(
    combined,
    format    = "latex",
    booktabs  = TRUE,
    linesep   = "",
    escape    = FALSE,
    col.names = col.names
  ) %>%
    kable_styling(
      latex_options = c("hold_position"),
      full_width    = FALSE,
      font_size     = 10
    ) %>%
    column_spec(1, width = "6em")

  cur <- 1
  for (i in seq_along(labels)) {
    tbl <- pack_rows(tbl, labels[i], cur, cur + emm_n[i] - 1)
    cur <- cur + emm_n[i]
  }
  tbl
}


# -------------------------------------------------------------------------
# Egg-count assay
# -------------------------------------------------------------------------
count_lm  <- lmer(Total   ~ Nuc * Treatment + Egg_Count + (1|Nuc:Mito), data = egg_count)
countF_lm <- lmer(F_total ~ Nuc * Treatment + Egg_Count + (1|Nuc:Mito), data = egg_count)
countM_lm <- lmer(M_total ~ Nuc * Treatment + Egg_Count + (1|Nuc:Mito), data = egg_count)

anova_table(
  models = list(count_lm, countF_lm, countM_lm),
  labels = c("Total", "Female", "Male")
)

# ANOVA for Egg_Count as a response — tests whether egg production itself
# varies by Mito, Nuc, and Treatment
countE_lm <- lmer(Egg_Count ~  Nuc * Treatment + (1|Nuc:Mito), data = egg_count)

anova_table(
  models = list(countE_lm),
  labels = c("Egg Count")
)

emm_table(
  models    = list(count_lm, countF_lm, countM_lm),
  labels    = c("Total", "Female", "Male"),
  spec      = ~ Treatment + Egg_Count,
  col.names = c("Treatment", "Egg Count", "EMM", "SE", "$df$", "Lower CI", "Upper CI")
)
