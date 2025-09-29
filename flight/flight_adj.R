############################################################
# Script: flight_adj.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Apply within-set and across-set corrections to 
#          raw flight assay data, generating adjusted 
#          phenotypes for downstream analysis and plotting.
############################################################

# -----------------------------------------
# Load required libraries
# -----------------------------------------
packages <- c("dplyr", "ggplot2")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

# -----------------------------------------
# Load flight data
# -----------------------------------------
flight = read.csv("flight/data/flight.csv", header = TRUE)

# -----------------------------------------
# Get within-set correction term
# -----------------------------------------
# For each Set × Sex × Treatment group:
#   - Keep only "yak" and "sm21" mito lines with Build A, or parentals
#   - Compute the mean phenotype (Y) as the within-set correction
phenos_correct = flight %>% 
  na.omit() %>%
  group_by(Set, Sex, Treatment) %>%
  filter((Mito %in% c("yak", "sm21") & Build == "A") | Build == "parental") %>%
  summarise(Set_adj = mean(Y), .groups = "drop")

# ---------------------------------------
# Get across-set correction term
# ---------------------------------------
# For each Treatment × Sex combination:
#   - Compute the mean phenotype (Y) as the across-set correction
pos_correct = flight %>% 
  na.omit() %>%
  group_by(Treatment, Sex) %>%
  summarise(Pos_adj = mean(Y), .groups = "drop")


# --------------------------------------------
# Apply correction terms
# --------------------------------------------
# Merge correction terms back into the full dataset
# Adjusted phenotype:
#   Y_adj = (Y - Set_adj) + Pos_adj
flight_adj = flight %>% 
  inner_join(phenos_correct, by = c("Set", "Sex", "Treatment")) %>%
  inner_join(pos_correct, by = c("Treatment", "Sex")) %>%
  mutate(Y_adj = (Y - Set_adj) + Pos_adj)

# --------------------------------------------
# Summarise adjusted phenotypes for plotting
# --------------------------------------------
# For each Mito × Nuc × Set × Sex × Treatment × Build:
#   - Compute the mean adjusted phenotype
plot_flight = flight_adj %>%
  group_by(Mito, Nuc, Set, Sex, Treatment, Build) %>%
  summarise(Y = mean(Y_adj), .groups = "drop")


# -------------------------------------------
# Save outputs
# -------------------------------------------
write.csv(flight_adj, "flight/dataflight_adj.csv", row.names = FALSE, quote = FALSE)

