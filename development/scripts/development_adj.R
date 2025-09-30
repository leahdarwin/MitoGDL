############################################################
# Script: development_adj.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Normalize pupation assay data by applying
#          within-set and across-set correction terms,
#          adjusting for larval density, and exporting
#          a cleaned dataset for downstream analyses.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("dplyr", "ggplot2")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------------------------------------
## Load development dataset
## ---------------------------------------------------------
dev <- read.csv("development/data/development.csv") %>%
  na.omit()  # remove missing values

## ---------------------------------------------------------
## Compute within-set correction terms
## - Restrict to reference genotypes (yak, sm21, or parental)
## - Calculate mean TSE for each Set × Treatment
## ---------------------------------------------------------
phenos_correct <- dev %>%
  group_by(Set, Treatment) %>%
  filter((Mito %in% c("yak", "sm21") & Build == "A") | Build == "parental") %>%
  summarise(Set_adj = mean(TSE), .groups = "drop")

## ---------------------------------------------------------
## Compute across-set correction terms
## - Mean TSE for each Treatment (pooled across sets)
## ---------------------------------------------------------
pos_correct <- dev %>%
  group_by(Treatment) %>%
  summarise(Pos_adj = mean(TSE), .groups = "drop")

## ---------------------------------------------------------
## Compute treatment-level total counts
## - Used to standardize larval density
## ---------------------------------------------------------
total_treat <- dev %>%
  group_by(Treatment) %>%
  summarise(Total_treat = mean(Total), .groups = "drop")

## ---------------------------------------------------------
## Apply correction terms
## - Adjust TSE: (TSE – Set_adj) + Pos_adj
## - Compute relative larval density
## ---------------------------------------------------------
dev_adj <- dev %>%
  inner_join(phenos_correct, by = c("Set", "Treatment")) %>%
  inner_join(pos_correct, by = "Treatment") %>%
  inner_join(total_treat, by = "Treatment") %>%
  mutate(
    Y_adj = (TSE - Set_adj) + Pos_adj,
    Larval_density = Total / Total_treat
  )

## ---------------------------------------------------------
## Export adjusted dataset
## ---------------------------------------------------------
write.csv(dev_adj, "development/data/development_adj.csv", row.names = FALSE, quote = FALSE)

