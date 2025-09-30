############################################################
# Script: weight_adj.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Process fly weight data by removing outliers,
#          applying within- and across-set corrections, and
#          generating an adjusted dataset for downstream 
#          analyses.
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
## Load raw weight data
## ---------------------------------------------------------
weight <- read.csv("weight/data/weight.csv", header = TRUE)

## ---------------------------------------------------------
## Count number of observations per group
## - Grouped by mitochondrial genotype, nuclear genotype,
##   treatment, and sex
## - Also compute mean per-fly weight per group
## ---------------------------------------------------------
counts <- weight %>%
  group_by(Mito, Nuc, Treatment, Sex) %>%
  summarise(count = n(), 
            Y = mean(weight_per_fly), .groups = "drop")

## ---------------------------------------------------------
## Remove outliers using z-scores within each sex
## - Calculate z-scores of per-fly weight
## - Keep observations within ±3 SD
## ---------------------------------------------------------
weight <- weight %>%
  group_by(Sex) %>%
  mutate(z = scale(weight_per_fly)) %>%
  filter(abs(z) <= 3) %>%
  ungroup()

## ---------------------------------------------------------
## Compute within-set correction terms
## - For each set, sex, and treatment
## - Use reference genotypes (yak/sm21 with Build A or parentals)
## - Correction = mean weight of reference group
## ---------------------------------------------------------
phenos_correct <- weight %>% 
  na.omit() %>%
  group_by(Set, Sex, Treatment) %>%
  filter((Mito %in% c("yak", "sm21") & Build == "A") | Build == "parental") %>%
  summarise(Set_adj = mean(weight_per_fly), .groups = "drop")

## ---------------------------------------------------------
## Compute across-set correction terms
## - For each treatment × sex combination
## - Correction = mean weight across sets
## ---------------------------------------------------------
pos_correct <- weight %>% 
  na.omit() %>%
  group_by(Treatment, Sex) %>%
  summarise(Pos_adj = mean(weight_per_fly), .groups = "drop")

## ---------------------------------------------------------
## Apply correction terms
## - Adjust Y by subtracting set-level correction and adding
##   across-set correction
## - Drops parental controls and set 0 during correction
## ---------------------------------------------------------
weight_adj <- weight %>% 
  inner_join(phenos_correct, join_by(Set, Sex, Treatment)) %>%
  inner_join(pos_correct, join_by(Treatment, Sex)) %>%
  mutate(Y_adj = (weight_per_fly - Set_adj) + Pos_adj)

## ---------------------------------------------------------
## Save adjusted dataset
## ---------------------------------------------------------
write.csv(weight_adj, "weight/data/weight_adj.csv", row.names = FALSE, quote = FALSE)
