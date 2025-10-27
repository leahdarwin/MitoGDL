############################################################
# Script: 
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: 
############################################################

# Load libraries ----------------------------------------------------------
packages <- c("dplyr", "ggplot2", "patchwork","lme4","emmeans",  "ggpubr", "ggtext")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

# -------------------------------------------------------------------
# Load datasets
# -------------------------------------------------------------------
weight <- read.csv("weight/data/weight_adj.csv") 
climb  <- read.csv("climb/data/climb_adj.csv") %>% filter(Build!="parental")
dev    <- read.csv("development/data/development_adj.csv")%>% filter(Build!="parental")

## Define consistent color palette for mitochondrial origins
color_palette <- c(
  "Beijing"   = "#1C448E",
  "Zimbabwe"  = "#52C2BA",
  "D.yakuba"  = "#FCAB10",
  "D.simulans"= "#ED1C24",
  "parent"    = "black"
)

plot_contrast = function(df){
  ggplot(df, aes(x = Mito, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = estimate-SE, ymax = estimate+SE), width = 0.1) +
    theme_classic() 
}

plot_ixn_mtnuc = function(df, label_mts){
  
  df = df %>%
    mutate(mitoOrig = case_when(
      grepl("B", Mito)      ~ "Beijing",
      grepl("Z", Mito)      ~ "Zimbabwe",
      grepl("siI", Mito)   ~ "D.simulans",
      grepl("yak", Mito)    ~ "D.yakuba",
      .default              = Mito
    ))
  
  labels_df = df %>%
    filter(Mito %in% label_mts & Nuc == "Ore")
  
  ggplot(df, aes(x = Nuc, y = emmean, group = Mito, color=mitoOrig)) +
   # geom_line(size = 1,alpha = ifelse(emm_df$Mito %in% label_mts, 1, 0.7)) +
    geom_line(size=1)+
    geom_point(size = 2) +
    theme_classic() +  
    scale_color_manual(values = setNames(color_palette[df$mitoOrig], df$mitoOrig)) + 
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))
    # + geom_text(data = labels_df,
    #   aes(label = Mito),
    #   vjust = -0.5,
    #   size = 3)
}


plot_ixn_mttreat = function(df, label_mts){
  
  df = df %>%
    mutate(mitoOrig = case_when(
      grepl("B", Mito)      ~ "Beijing",
      grepl("Z", Mito)      ~ "Zimbabwe",
      grepl("siI", Mito)   ~ "D.simulans",
      grepl("yak", Mito)    ~ "D.yakuba",
      .default              = Mito
    ))
  
  labels_df = df %>%
    filter(Mito %in% label_mts & Treatment == "Rotenone")
  
  ggplot(df, aes(x = Treatment, y = emmean, group = Mito, color=mitoOrig)) +
  #  geom_line(size = 1,alpha = ifelse(emm_df$Mito %in% label_mts, 1, 0.7)) +
    geom_line(size = 1)+
    geom_point(size = 2) +
    theme_classic() + 
    scale_color_manual(values = setNames(color_palette[df$mitoOrig], df$mitoOrig)) + 
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)))
    # + geom_text(data = labels_df,
    #           aes(label = Mito),
    #           vjust = -0.5,
    #           size = 3)
}

get_dfs = function(emm){
  
  contrast = contrast(emm, method = "pairwise")
  contrast_df = as.data.frame(contrast)
  emm_df = as.data.frame(emm) 
  
  return(list(contrast_df, emm_df))
}

##Filter weight
weightM = weight %>% filter(Sex=="M", Build!="parental") %>% mutate(Y_adj = Y_adj*1000)

##Define weight model
weightM_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build),
  data = weightM
)

##Estimate marginal means
emm_mito_weight = emmeans(weightM_lm, ~ Nuc | Mito)
weight_dfs = get_dfs(emm_mito_weight)

weight_contrast_plot = plot_contrast(weight_dfs[[1]])
weight_contrast_plot
weight_ixn_plot = plot_ixn_mtnuc(weight_dfs[[2]], c(weight_dfs[[1]]$Mito[which.max(weight_dfs[[1]]$estimate)], weight_dfs[[1]]$Mito[which.min(weight_dfs[[1]]$estimate)])) +
  labs( x = "Nuclear Background",
        y = "EM mean weight (mg)",
        title = "") + theme(legend.position = "none") + 
  theme(
    axis.title.y = element_textbox_simple(
      size = 16, lineheight = 1, halign = 0.5,
      fill = "#fbb4ae", padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 5, 0), 
      orientation = "left-rotated"
    ))
weight_ixn_plot 

dev_lm <- lmerTest::lmer(
  Y_adj ~ Mito * Nuc * Treatment + Larval_density + (1 | Mito:Nuc:Build),
  data = dev
)

##Estimate marginal means
emm_mito_dev = emmeans(dev_lm, ~ Treatment | Mito)
dev_dfs = get_dfs(emm_mito_dev)

dev_contrast_plot = plot_contrast(dev_dfs[[1]])
dev_contrast_plot
dev_ixn_plot = plot_ixn_mttreat(dev_dfs[[2]], c(dev_dfs[[1]]$Mito[which.max(dev_dfs[[1]]$estimate)], dev_dfs[[1]]$Mito[which.min(dev_dfs[[1]]$estimate)])) +
  labs( x = "Environment",
        y = "EM mean development time (days)",
        title = "") + theme(legend.position = "none") +
  theme(
    axis.title.y = element_textbox_simple(
      size = 16, lineheight = 1, halign = 0.5,
      fill = "#decbe4", padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 5, 0), 
      orientation = "left-rotated"
    ))
dev_ixn_plot


climbM = climb %>% filter(Sex=="M")
climbM_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) + 
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbM
)

emm_mito_climbM = emmeans(climbM_lm, ~ Treatment | Mito)
climbM_dfs = get_dfs(emm_mito_climbM)

climbM_contrast_plot = plot_contrast(climbM_dfs[[1]])
climbM_contrast_plot
climbM_ixn_plot = plot_ixn_mttreat(climbM_dfs[[2]], c(climbM_dfs[[1]]$Mito[which.max(climbM_dfs[[1]]$estimate)], climbM_dfs[[1]]$Mito[which.min(climbM_dfs[[1]]$estimate)])) +
  labs( x = "Environment",
        y = "",
        title = "") + theme(legend.position = "none")
climbM_ixn_plot

climbF = climb %>% filter(Sex=="F")
climbF_lm <- lmer(
  Y_adj ~ Mito * Nuc * Treatment + (1 | Mito:Nuc:Build) + 
    (1 | Mito:Nuc:Build:Treatment:Vial),
  data = climbF
)

emm_mito_climbF = emmeans(climbF_lm, ~ Treatment | Mito)
climbF_dfs = get_dfs(emm_mito_climbF)

climbF_contrast_plot = plot_contrast(climbF_dfs[[1]])
climbF_contrast_plot
climbF_ixn_plot = plot_ixn_mttreat(climbF_dfs[[2]], c(climbF_dfs[[1]]$Mito[which.max(climbF_dfs[[1]]$estimate)], climbF_dfs[[1]]$Mito[which.min(climbF_dfs[[1]]$estimate)])) +
  labs( x = "Environment",
        y = "EM mean climbing speed (cm/s)",
        title = "") + theme(legend.position = "none") +
  theme(
    axis.title.y = element_textbox_simple(
      size = 16, lineheight = 1, halign = 0.5,
      fill = "#b3cde3", padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 5, 0), 
      orientation = "left-rotated"
    ))
  
climbF_ixn_plot

weight_ixn_plot  | dev_ixn_plot | climbF_ixn_plot | climbM_ixn_plot
##save as 11x4.5
