###########################################################
# Script: mainEffects_plot.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Plots emmeans and adjusted averages for significant main effects. 
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
weight = read.csv("weight/data/weight_adj.csv") %>% filter(Build != "parental")
climb  = read.csv("climb/data/climb_adj.csv") %>% filter(Build != "parental")
dev    = read.csv("development/data/development_adj.csv") %>% filter(Build != "parental")
flight = read.csv("flight/data/flight_adj.csv") %>% filter(Build != "parental")

# ===================================================================
# Helper plotting functions
# ===================================================================

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

## ---------------------------------------------------------------
## Plot treatment means (raw data + model-adjusted EMMs)
## ---------------------------------------------------------------
plot_treatment_means = function(df, lm) {
  
  # Get estimated marginal means for Treatment
  emm_df = as.data.frame(emmeans(lm, ~ Treatment))
  
  # Summarize observed means (for raw data visualization)
  plot_df = df %>%
    group_by(Mito, Nuc, Treatment, Build) %>%
    summarise(Y_avg = mean(Y_adj), .groups = "drop") %>%
    na.omit()
  
  # Plot: raw jittered points + model means with SE
  ggplot(plot_df, aes(x = Treatment, y = Y_avg, color = Treatment)) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4),
      alpha = 0.4, size = 1.5
    ) +
    geom_point(
      data = emm_df, aes(y = emmean),
      position = position_dodge(width = 0.4),
      size = 3, shape = 21, fill = "white"
    ) +
    theme_classic() +
    scale_color_manual(values = c("Control" = "black", "Rotenone" = "#723a83")) +
    theme(legend.position = "none") +
    labs(x = "")
}


## ---------------------------------------------------------------
## Plot nuclear background means (raw data + model-adjusted EMMs)
## ---------------------------------------------------------------
plot_nuc_means = function(df, lm) {
  
  # Get estimated marginal means for Nuclear background
  emm_df = as.data.frame(emmeans(lm, ~ Nuc))
  
  # Summarize observed means (for raw data visualization)
  plot_df = df %>%
    group_by(Mito, Nuc, Treatment, Build) %>%
    summarise(Y_avg = mean(Y_adj), .groups = "drop")%>%
    na.omit()
  
  # Plot: raw jittered points + model means with SE
  ggplot(plot_df, aes(x = Nuc, y = Y_avg, color = Nuc)) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4),
      alpha = 0.4, size = 1.5
    ) +
    geom_point(
      data = emm_df, aes(y = emmean),
      position = position_dodge(width = 0.4),
      size = 3, shape = 21, fill = "white"
    ) +
    theme_classic() +
    scale_color_manual(values = c("Ore" = "#67a9cf", "375" = "#ef8a62")) +
    theme(legend.position = "none") +
    labs(x = "")
}

## ---------------------------------------------------------------
## Plot mito means (raw data + model-adjusted EMMs)
## ---------------------------------------------------------------
plot_mito_means = function(df, lm) {
  
  # Get estimated marginal means for Treatment
  emm_df = as.data.frame(emmeans(lm, ~ Mito)) %>% 
    mutate(mitoOrig = map_mito_origin(Mito)) %>%
    arrange(emmean) %>%
    mutate(Mito = factor(Mito, levels = unique(Mito)))
  
  # Summarize observed means (for raw data visualization)
  plot_df = df %>%
    group_by(Mito, Nuc, Treatment, Build) %>%
    summarise(Y_avg = mean(Y_adj), .groups = "drop") %>% 
    mutate(mitoOrig = map_mito_origin(Mito))  %>%
    mutate(Mito = factor(Mito, levels = levels(emm_df$Mito))) %>%
    na.omit()
  
  # Plot: raw jittered points + model means with SE
  ggplot(plot_df, aes(x = Mito, y = Y_avg, color = mitoOrig)) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4),
      alpha = 0.4, size = 1.5
    ) +
    geom_point(
      data = emm_df, aes(y = emmean),
      position = position_dodge(width = 0.4),
      size = 3, shape = 21, fill = "white"
    ) +
    theme_classic() +
    scale_color_manual(values = color_palette) +
    theme(legend.position = "none") +
    labs(x = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
}

# ===================================================================
# Models and plots for each trait
# ===================================================================

## --------------------------
## Weight (Females)
## --------------------------
weightF = weight %>%
  filter(Sex == "F") %>%
  mutate(Y_adj = Y_adj * 1000)  # convert to mg

weightF_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = weightF
)

weightF_plot_nuc = plot_nuc_means(weightF, weightF_lm) +
  labs(y = "Weight (mg)")

weightF_plot_mt = plot_mito_means(weightF, weightF_lm) +
  labs(y = "Weight (mg)")

## --------------------------
## Weight (Males)
## --------------------------
weightM = weight %>%
  filter(Sex == "M") %>%
  mutate(Y_adj = Y_adj * 1000)  # convert to mg

weightM_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = weightM
)

weightM_plot_mt = plot_mito_means(weightM, weightM_lm) +
  labs(y = "Weight (mg)")

## --------------------------
## Climbing (Males)
## --------------------------
climbM = climb %>% filter(Sex == "M")

climbM_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment +
    (1 | Mito:Nuc:Build) +
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbM
)

climbM_plot_treat = plot_treatment_means(climbM, climbM_lm) +
  labs(y = "Climbing speed (cm/s)")

climbM_plot_nuc = plot_nuc_means(climbM, climbM_lm) +
  labs(y = "Climbing speed (cm/s)")


## --------------------------
## Climbing (Females)
## --------------------------
climbF = climb %>% filter(Sex == "F")

climbF_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment +
    (1 | Mito:Nuc:Build) +
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbF
)

climbF_plot_treat = plot_treatment_means(climbF, climbF_lm) +
  labs(y = "Climbing \nspeed (cm/s)")

climbF_plot_nuc = plot_nuc_means(climbF, climbF_lm) +
  labs(y = "Climbing \nspeed (cm/s)")

climbF_plot_mt = plot_mito_means(climbF, climbF_lm) +
  labs(y = "Climbing \nspeed (cm/s)") 

## --------------------------
## Flight (Males)
## --------------------------
flightM = flight %>% filter(Sex == "M")

flightM_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = flightM
)

flightM_plot_treat = plot_treatment_means(flightM, flightM_lm) +
  labs(y = "Flight landing \nheight (m)")

flightM_plot_nuc = plot_nuc_means(flightM, flightM_lm) +
  labs(y = "Flight landing \nheight (m)")


## --------------------------
## Flight (Females)
## --------------------------
flightF = flight %>% filter(Sex == "F")

flightF_lm = lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = flightF
)

flightF_plot_nuc = plot_nuc_means(flightF, flightF_lm) +
  labs(y = "Flight landing \nheight (m)")


## --------------------------
## Development
## --------------------------
dev_lm = lmerTest::lmer(
  Y_adj ~ Mito * Nuc * Treatment + Larval_density + (1 | Mito:Nuc:Build),
  data = dev
)

dev_plot_treat = plot_treatment_means(dev, dev_lm) +
  labs(y = "Development \ntime (days)")

dev_plot_nuc = plot_nuc_means(dev, dev_lm) +
  labs(y = "Development \ntime (days)")


# ===================================================================
# Combine plots into panels
# ===================================================================

## ---------------------------------------------------------------
## Treatment effects
## ---------------------------------------------------------------
treat = (
  flightM_plot_treat + 
    dev_plot_treat + 
    climbF_plot_treat + ylim(c(0.3, 4.8)) +
    climbM_plot_treat + ylim(c(0.3, 4.8)) +
    labs(y = "")
) & 
  scale_x_discrete(labels = c("Control" = "C", "Rotenone" = "R"))

treat
## ---------------------------------------------------------------
## Nuclear background effects
## ---------------------------------------------------------------
nuc = (
  flightF_plot_nuc + ylim(c(0.55, 0.91)) +
    flightM_plot_nuc + ylim(c(0.55, 0.91)) + labs(y = "") +
    dev_plot_nuc +
    climbF_plot_nuc + ylim(c(0.4, 4.8)) +
    climbM_plot_nuc + ylim(c(0.4, 4.8)) + labs(y = "") +
    weightF_plot_nuc
)

nuc
# ===================================================================
# Combine all panels
# ===================================================================
nuc_treat = (nuc | treat) + 
  plot_layout(widths = c(1, 0.66), guides = "collect")
mito = weightF_plot_mt / weightM_plot_mt / climbF_plot_mt

final = nuc_treat / mito +
  plot_layout(heights = c(0.6,1))

ggsave("figures/main_figs/mainEffects_fig2.pdf", final, width = 6.5, height = 7.5)

