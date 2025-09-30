############################################################
# Script: phenotype_pca.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Summarize fly phenotype data (weight, climb,
#          flight, development) by sex, treatment, and genotype,
#          merge into a combined dataset, and run PCA to 
#          visualize multivariate trait patterns.
############################################################

# -------------------------------------------------------------------
# Load libraries
# -------------------------------------------------------------------
packages <- c("dplyr", "ggplot2", "ggfortify", "ggsci", "grid", "gridExtra")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p, dependencies = TRUE)
}
lapply(packages, library, character.only = TRUE)

# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
weight <- read.csv("weight/data/weight.csv")
climb  <- read.csv("climb/data/climb.csv")
dev    <- read.csv("development/data/development.csv")
flight <- read.csv("flight/data/flight.csv")


# -------------------------------------------------------------------
# Summarize by sex and trait
# -------------------------------------------------------------------
# Weight
weight_F <- weight %>%
  na.omit() %>%
  filter(Sex == "F") %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(weight_F = mean(weight_per_fly))

weight_M <- weight %>%
  na.omit() %>%
  filter(Sex == "M") %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(weight_M = mean(weight_per_fly))

# Climb
climb_F <- climb %>%
  na.omit() %>%
  filter(Sex == "F") %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(climb_F = mean(Slope))

climb_M <- climb %>%
  na.omit() %>%
  filter(Sex == "M") %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(climb_M = mean(Slope))

# Flight
flight_F <- flight %>%
  na.omit() %>%
  filter(Sex == "F") %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(flight_F = mean(Y))

flight_M <- flight %>%
  na.omit() %>%
  filter(Sex == "M") %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(flight_M = mean(Y))

# Development (no sex split)
dev_MF <- dev %>%
  na.omit() %>%
  group_by(Mito, Nuc, Treatment) %>%
  summarise(dev = mean(TSE))

# -------------------------------------------------------------------
# Merge datasets
# -------------------------------------------------------------------
merged <- full_join(dev_MF,
                    full_join(weight_M,
                              full_join(climb_F,
                                        full_join(climb_M,
                                                  full_join(flight_F,
                                                            full_join(flight_M, weight_F)))))) %>%
  na.omit()

# -------------------------------------------------------------------
# Principal Component Analysis (PCA)
# -------------------------------------------------------------------
# Overall PCA
pca <- prcomp(merged[, 4:10], scale. = TRUE)

# -------------------------------------------------------------------
# Plot PCA
# -------------------------------------------------------------------
p = autoplot(pca, data = merged, color = "Nuc", x = 1, y = 2) +
  scale_color_manual(values = c("Ore" = "#67a9cf", "375" = "#ef8a62")) +
  theme_light() +
  labs(color = "Nuclear Genotype") + 
  theme(legend.position = "top")

# Save to PDF
ggsave("figures/pca_figS2.pdf", p, width = 8, height = 5)

