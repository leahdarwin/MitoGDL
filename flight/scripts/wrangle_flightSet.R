############################################################
## Script: wrangle_flightSet.R
## Author: Leah Darwin
## Date: 2025-09-29
## Description:
##   - Reads in raw flight assay CSV files
##   - Cleans and harmonizes column formats
##   - Assigns Set, Build, and Treatment values
##   - Joins genotype metadata from stock_genotype.csv
##   - Outputs final cleaned dataset: flight.csv
############################################################


## ---------------------------
## Load required packages
## ---------------------------
packages <- c("dplyr", "stringr")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------
## Input files
## ---------------------------

# Find all CSVs beginning with "flightSet"
files <- list.files(path = "./flight/data/", pattern = "^flightSet.*csv$", full.names = TRUE)
dfs   <- lapply(files, read.csv)

# Genotype metadata (drop precomputed MitoNuc column, recode Build)
stock_geno <- read.csv("stock_genotype.csv") %>% 
  select(-MitoNuc) %>%
  mutate(Build = case_when(
    grepl("A", Stock) ~ "A",
    grepl("B", Stock) ~ "B",
    .default = "parental"
  ))

## ---------------------------
## Define experiment metadata
## ---------------------------
# Each flight set corresponds to a Set, Build, and Treatment
sets   <- c(1, 1, 2, 2, 3, 3, 4, 5, 6, 7)
builds <- c(NA, NA, "A", "B", "A", "B", NA, NA, NA, NA)
treats <- c("Control", "Rotenone", rep(NA, 8))


## ---------------------------
## General edits applied to all dataframes
## ---------------------------
dfs <- lapply(seq_along(dfs), function(i) {
  dfs[[i]] %>%
    select(Mito, Nuc, Sex, Y) %>%                  # keep relevant columns
    mutate(Set = sets[i],
           Build = builds[i],
           Treatment = treats[i]) %>%
    filter(!grepl("test", Nuc)) %>%                # drop test rows
    filter(!grepl("w1118", Nuc)) %>%               # drop unwanted control
    filter(Y > 0)                                  # keep only positive phenotypes
})


## ---------------------------
## Set-specific cleaning rules
## ---------------------------

# Sets 1–2: parental vs A/B builds
for (i in c(1, 2)) {
  dfs[[i]] <- dfs[[i]] %>%
    mutate(Build = case_when(
      Nuc %in% c("OreR", "375")   ~ "parental",
      Nuc %in% c("OreB", "375B") ~ "B",
      Nuc %in% c("OreA", "375A") ~ "A"
    ),
    Nuc = case_when(
      grepl("Ore", Nuc) ~ "Ore",
      grepl("375", Nuc) ~ "375"
    ))
}

# Sets 3–6: treatment encoded in Nuc suffix
for (i in c(3:6)) {
  dfs[[i]] <- dfs[[i]] %>%
    mutate(Treatment = case_when(
      grepl("R$", Nuc) ~ "Rotenone",
      grepl("C$", Nuc) ~ "Control"
    ),
    Nuc = case_when(
      grepl("Ore", Nuc) ~ "Ore",
      grepl("375", Nuc) ~ "375"
    ))
}

# Sets 7–10: treatment + build encoded in Nuc suffix
for (i in c(7:10)) {
  dfs[[i]] <- dfs[[i]] %>%
    mutate(Treatment = case_when(
      grepl("R$", Nuc) ~ "Rotenone",
      grepl("C$", Nuc) ~ "Control"
    ),
    Build = case_when(
      grepl("A", Nuc) ~ "A",
      grepl("B", Nuc) ~ "B",
      .default        = "parental"
    ),
    Nuc = case_when(
      grepl("Ore", Nuc) ~ "Ore",
      grepl("375", Nuc) ~ "375"
    ))
}

## ---------------------------
## Combine and finalize dataset
## ---------------------------
df <- bind_rows(dfs) %>%
  # Standardize mitochondrial haplotype names
  mutate(Mito = gsub("A", "", Mito),
         Mito = gsub("B$", "", Mito),
         Mito = gsub("R$", "", Mito),
         Mito = case_when(
           Mito %in% c("parental", "Parental") ~ Nuc,
           .default                            = Mito
         )) %>%
  # Add genotype metadata
  left_join(stock_geno, join_by(Mito, Nuc, Build), relationship = "many-to-one") %>%
  filter(Mito != "Bei") %>%
  # Ensure parental Build is properly coded
  mutate(Build = case_when(
    Mito == Nuc ~ "parental",
    .default    = Build
  ))

## ---------------------------
## Write output
## ---------------------------
write.csv(df, "flight/data/flight.csv", row.names = FALSE, quote = FALSE)
