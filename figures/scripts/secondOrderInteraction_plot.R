###########################################################
# Script: secondOrderInteraction_plot.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Plots interaction plots for significant second order effects using estimated marginal means. 
############################################################

# Load libraries ----------------------------------------------------------
## -------------------------------------------------------------------
## Setup: load & install required packages
## -------------------------------------------------------------------
packages = c(
  "dplyr", "ggplot2", "patchwork",
  "lme4", "lmerTest", "emmeans",
  "ggpubr", "ggtext"
)

new_pkgs = setdiff(packages, rownames(installed.packages()))
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, dependencies = TRUE)
}
invisible(lapply(packages, library, character.only = TRUE))

## -------------------------------------------------------------------
## Load datasets
## -------------------------------------------------------------------
weight = read.csv("weight/data/weight_adj.csv")
climb  = read.csv("climb/data/climb_adj.csv") %>% filter(Build != "parental")
dev    = read.csv("development/data/development_adj.csv") %>% filter(Build != "parental")

## -------------------------------------------------------------------
## Define consistent color palette
## -------------------------------------------------------------------
color_palette = c(
  "Beijing"    = "#1C448E",
  "Zimbabwe"   = "#52C2BA",
  "D.yakuba"   = "#FCAB10",
  "D.simulans" = "#ED1C24",
  "parent"     = "black"
)

## -------------------------------------------------------------------
## Helper functions
## -------------------------------------------------------------------

# Map mitochondrial codes to full origin names
map_mito_origin = function(mito_vector) {
  dplyr::case_when(
    grepl("B", mito_vector)      ~ "Beijing",
    grepl("Z", mito_vector)      ~ "Zimbabwe",
    grepl("siI", mito_vector)    ~ "D.simulans",
    grepl("yak", mito_vector)    ~ "D.yakuba",
    TRUE                         ~ mito_vector
  )
}

# Extract both emmeans and contrasts data frames
get_dfs = function(emm_obj) {
  contr = contrast(emm_obj, method = "pairwise")
  list(
    contrasts = as.data.frame(contr),
    emmeans   = as.data.frame(emm_obj)
  )
}

# Simple contrast plot (mean ± SE by mito)
plot_contrast = function(df) {
  ggplot(df, aes(x = Mito, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.1) +
    theme_classic() +
    xlab("Mito") +
    ylab("Estimate")
}

# Mito × Nuc interaction plot
plot_ixn_mtnuc = function(df, palette = color_palette) {
  df2 = df %>% mutate(mitoOrig = map_mito_origin(Mito))
  ggplot(df2, aes(x = Nuc, y = emmean, group = Mito, color = mitoOrig)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_classic() +
    scale_color_manual(values = palette) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    labs(x = "Nuclear Background", y = NULL)
}

# Mito × Treatment interaction plot
plot_ixn_mttreat = function(df, palette = color_palette) {
  df2 = df %>% mutate(mitoOrig = map_mito_origin(Mito))
  ggplot(df2, aes(x = Treatment, y = emmean, group = Mito, color = mitoOrig)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_classic() +
    scale_color_manual(values = palette) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    labs(x = "Environment", y = NULL)
}

## -------------------------------------------------------------------
## ANALYSES & PLOTS
## -------------------------------------------------------------------

### --- Weight (Males) ---
weightM = weight %>%
  filter(Sex == "M", Build != "parental") %>%
  mutate(Y_adj = Y_adj * 1000)

weightM_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = weightM
)

emm_mito_weight = emmeans(weightM_lm, ~ Nuc | Mito)
weight_dfs = get_dfs(emm_mito_weight)

weight_contrast_plot = plot_contrast(weight_dfs$contrasts)
weight_ixn_plot = plot_ixn_mtnuc(weight_dfs$emmeans) +
  labs(y = "EM mean weight (mg)") +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_textbox_simple(
      size = 16, lineheight = 1, halign = 0.5,
      fill = "#fbb4ae", padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 5, 0),
      orientation = "left-rotated"
    )
  )

### --- Development ---
dev_lm = lmerTest::lmer(
  Y_adj ~ Mito * Nuc * Treatment + Larval_density + (1 | Mito:Nuc:Build),
  data = dev
)

emm_mito_dev = emmeans(dev_lm, ~ Treatment | Mito)
dev_dfs = get_dfs(emm_mito_dev)

dev_contrast_plot = plot_contrast(dev_dfs$contrasts)
dev_ixn_plot = plot_ixn_mttreat(dev_dfs$emmeans) +
  labs(y = "EM mean development time (days)") +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_textbox_simple(
      size = 16, lineheight = 1, halign = 0.5,
      fill = "#decbe4", padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 5, 0),
      orientation = "left-rotated"
    )
  )

### --- Climbing (Males) ---
climbM = climb %>% filter(Sex == "M")
climbM_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) +
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbM
)

emm_mito_climbM = emmeans(climbM_lm, ~ Treatment | Mito)
climbM_dfs = get_dfs(emm_mito_climbM)

climbM_contrast_plot = plot_contrast(climbM_dfs$contrasts)
climbM_ixn_plot = plot_ixn_mttreat(climbM_dfs$emmeans) +
  labs(y = NULL) +
  theme(legend.position = "none")

### --- Climbing (Females) ---
climbF = climb %>% filter(Sex == "F")
climbF_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) +
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbF
)

emm_mito_climbF = emmeans(climbF_lm, ~ Treatment | Mito)
climbF_dfs = get_dfs(emm_mito_climbF)

climbF_contrast_plot = plot_contrast(climbF_dfs$contrasts)
climbF_ixn_plot = plot_ixn_mttreat(climbF_dfs$emmeans) +
  labs(y = "EM mean climbing speed (cm/s)") +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_textbox_simple(
      size = 16, lineheight = 1, halign = 0.5,
      fill = "#b3cde3", padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 5, 0),
      orientation = "left-rotated"
    )
  )

## -------------------------------------------------------------------
## Combine and save plots
## -------------------------------------------------------------------
combined_plot = weight_ixn_plot | dev_ixn_plot | climbF_ixn_plot | climbM_ixn_plot
ggsave("figures/main_figs/secondOrderInteraction_fig3.pdf", combined_plot, width = 11, height = 4.5)



