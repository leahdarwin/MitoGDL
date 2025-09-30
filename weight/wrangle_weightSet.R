############################################################
# Script: wrangle_weightSet.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Process fly body weight data by merging multiple
#          sets, correcting annotations, computing per-fly
#          weight, and standardizing genotype labels.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("dplyr")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------------------------------------
## Load stock genotype reference
## ---------------------------------------------------------
stock_geno <- read.csv("stock_genotype.csv")

## ---------------------------------------------------------
## Load raw weight set CSVs
## - Matches files beginning with "weightSet"
## - Combine into a single dataframe
## ---------------------------------------------------------
files <- list.files(path = "weight/data/", pattern = "^weightSet.*csv$", full.names = TRUE)
dfs <- lapply(files, read.csv)
weight_df <- bind_rows(dfs)

## ---------------------------------------------------------
## Merge with genotype reference
## - Keep only rows with matching entries
## - Drop rows with missing data
## ---------------------------------------------------------
weight_df <- merge(weight_df, stock_geno, all = TRUE)
weight_df <- na.omit(weight_df)

## ---------------------------------------------------------
## Compute per-fly weight
## - Weight_flies = total vial weight
## - Weight       = empty vial weight
## - Inds         = number of flies
## - Filter out non-positive values
## ---------------------------------------------------------
weight_df$weight_per_fly <- (weight_df$Weight_flies - weight_df$Weight) / weight_df$Inds
weight_df <- subset(weight_df, weight_per_fly > 0)

## ---------------------------------------------------------
## Standardize Treatment names
## - Fix spacing in "Rotenone"
## ---------------------------------------------------------
weight_df$Treatment <- gsub("Rotenone ", "Rotenone", weight_df$Treatment)

## ---------------------------------------------------------
## Assign build categories
## - A, B, or parental (default)
## ---------------------------------------------------------
weight_df$Build <- ifelse(grepl("A", weight_df$Stock), "A",
                          ifelse(grepl("B", weight_df$Stock), "B",
                                 "parental"))

## ---------------------------------------------------------
## Standardize mitochondrial names
## - Capitalize "Zim53"
## ---------------------------------------------------------
weight_df <- weight_df %>%
  mutate(Mito = case_when(
    Mito == "zim53" ~ "Zim53",
    .default = Mito
  ))

## ---------------------------------------------------------
## Save cleaned dataset
## ---------------------------------------------------------
write.csv(weight_df, file = "weight/data/weight.csv", row.names = FALSE, quote = FALSE)
