############################################################
# Script: set9Corr_plot.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Compare split vs. unified experiments across
#          climbing, weight, and development phenotypes.
#          Produces rank–rank correlation plots and stats.
############################################################

# -------------------------------------------------------------------
# Load libraries
# -------------------------------------------------------------------
packages <- c("dplyr", "ggplot2", "patchwork")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
climb  <- read.csv("climb/data/climb_adj.csv") %>% rename(Y = Slope)
dev    <- read.csv("development/data/development_adj.csv") %>% 
  mutate(Sex = "MF") %>% rename(Y = TSE)
weight <- read.csv("weight/data/weight_adj.csv") %>% rename(Y = weight_per_fly)

# Stocks included in Set 9 (special unified experiment)
stocks <- unique(climb$Stock[climb$Set == 9])

# Colors for sex-specific plots
sex_colors <- c("F" = "lightgrey", "M" = "gray34")

# -------------------------------------------------------------------
# Function: Prepare merged dataset for split vs. unified experiments
# -------------------------------------------------------------------
get_merged_df <- function(df) {
  
  # Split experiment (exclude Set 9)
  not_set9 <- df %>%
    filter(Stock %in% stocks, Set != 9) %>%
    group_by(Sex, Mito, Nuc, Treatment) %>%
    summarise(Y = mean(Y_adj), .groups = "drop")
  
  # Unified experiment (Set 9 only)
  set9 <- df %>%
    filter(Set == 9) %>%
    group_by(Sex, Mito, Nuc, Treatment) %>%
    summarise(Y = mean(Y), .groups = "drop")
  
  # Merge and compute comparison metrics
  merged <- full_join(set9, not_set9, by = join_by(Treatment, Mito, Nuc, Sex)) %>%
    na.omit() %>%
    group_by(Sex) %>%
    mutate(
      avg    = (Y.x + Y.y) / 2,   # Average across experiments
      diff   = Y.y - Y.x,         # Difference between experiments
      rank.x = rank(-Y.x),        # Rank in unified experiment
      rank.y = rank(-Y.y)         # Rank in split experiment
    ) %>%
    ungroup()
  
  return(merged)
}

# -------------------------------------------------------------------
# Function: Compute correlations and return a ggplot
# -------------------------------------------------------------------
get_corr <- function(merged, trait) {
  
  if (trait %in% c("Weight", "Climbing")) {
    # Female-specific correlation
    mergedF <- merged %>% filter(Sex == "F")
    message(trait, " (females):")
    print(cor.test(mergedF$Y.x, mergedF$Y.y, method = "spearman"))
    
    # Male-specific correlation
    mergedM <- merged %>% filter(Sex == "M")
    message(trait, " (males):")
    print(cor.test(mergedM$Y.x, mergedM$Y.y, method = "spearman"))
    
    # Plot with sex colors
    p <- ggplot(merged, aes(x = rank.x, y = rank.y, color = Sex)) +
      geom_point(size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = sex_colors) +
      labs(
        title = trait,
        x = "Rank in Split Experiment",
        y = "Rank in Unified Experiment"
      ) +
      theme_minimal(base_size = 14)
    
  } else {
    # Sex not relevant (development)
    message(trait, ":")
    print(cor.test(merged$Y.x, merged$Y.y, method = "spearman"))
    
    # Plot without sex colors
    p <- ggplot(merged, aes(x = rank.x, y = rank.y)) +
      geom_point(size = 3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = trait,
        x = "Rank in Split Experiment",
        y = "Rank in Unified Experiment"
      ) +
      theme_minimal(base_size = 14)
  }
  
  return(p)
}

# -------------------------------------------------------------------
# Run functions for each trait
# -------------------------------------------------------------------
merged_climb  <- get_merged_df(climb)
merged_weight <- get_merged_df(weight)
merged_dev    <- get_merged_df(dev)

p1 <- get_corr(merged_climb,  "Climbing")
p2 <- get_corr(merged_weight, "Weight")
p3 <- get_corr(merged_dev,    "Development")

# Combine into a single figure
combined_plot <- p3 + p1 + theme(legend.position = "none") + p2

# Save to PDF
ggsave("figures/set9_corr_figS3.pdf", combined_plot, width = 11, height = 4)

