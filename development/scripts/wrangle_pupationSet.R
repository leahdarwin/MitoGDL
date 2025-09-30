############################################################
# Script: wrangle_pupationSet.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Process pupation assay data by parsing raw CSVs,
#          extracting genotype and metadata annotations,
#          computing time-since-egg (TSE),
#          and merging with stock genotypes to produce a
#          clean, standardized dataset for analysis.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("car", "ggplot2", "dplyr", "stringr")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------------------------------------
## Load stock genotype metadata
## ---------------------------------------------------------
stock_geno <- read.csv("stock_genotype.csv")

## ---------------------------------------------------------
## Load pupation data (all files matching "pupationSet*.csv")
## ---------------------------------------------------------
files <- list.files(path = "development/data", pattern = "^pupationSet.*csv$",full.names = TRUE)
dfs <- lapply(files, read.csv)

## ---------------------------------------------------------
## Annotate set numbers to each dataset
## Note: `sets` defines the numeric Set labels
## ---------------------------------------------------------
sets <- c(2:9, 9)
dfs <- lapply(seq_along(dfs), function(i) {
  dfs[[i]] %>%
    mutate(Set = sets[i])
})

## ---------------------------------------------------------
## Define weights for each developmental day (6.5 to 12.5)
## ---------------------------------------------------------
day_weights <- c(6.5:12.5)


## ---------------------------------------------------------
## Wrangling function
## - Computes TSE as weighted mean of emergence counts
## - Classifies Build (A, B, or parental) from Stock ID
## - Computes total count across days
## - Joins with stock genotype annotations
## ---------------------------------------------------------
wrangle <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      ## Time to successful emergence (weighted average)
      TSE = sum(c_across(3:9) * day_weights) / sum(c_across(3:9)),
      
      ## Annotate build (A, B, parental) from Stock ID
      build = case_when(
        str_detect(Stock, "A") ~ "A",
        str_detect(Stock, "B") ~ "B",
        TRUE ~ "parental"
      ),
      
      ## Total count across developmental days
      total = sum(c_across(3:9))
    ) %>%
    ungroup() %>%
    select(-c(3:9)) %>%              # drop daily count columns
    inner_join(stock_geno, by = "Stock")  # add genotype annotations
}

## ---------------------------------------------------------
## Apply wrangling function to all sets and merge
## ---------------------------------------------------------
wrangled_dfs <- lapply(dfs, wrangle)
all_sets <- bind_rows(wrangled_dfs)

## ---------------------------------------------------------
## Standardize column names and annotations
## ---------------------------------------------------------
colnames(all_sets) <- c(
  "Stock", "Vial", "Treatment", "Set",
  "TSE", "Build", "Total", "MitoNuc", "Mito", "Nuc"
)

## Ensure Mito naming consistency (e.g., Zim53 capitalization)
all_sets <- all_sets %>%
  mutate(Mito = case_when(
    Mito == "zim53" ~ "Zim53",
    .default = Mito
  ))

## ---------------------------------------------------------
## Write final cleaned dataset
## ---------------------------------------------------------
write.csv(all_sets, file = "development/data/development.csv", row.names = FALSE, quote = FALSE)


