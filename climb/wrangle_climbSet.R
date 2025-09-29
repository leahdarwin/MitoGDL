############################################################
# Script: wrangle_climbSet.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Process climbing assay data by cleaning, merging,
#          and standardizing genotype annotations to create
#          a reproducible dataset for downstream analyses.
############################################################

required_packages <- c("dplyr", "stringr")

# Install any missing packages
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))


#----------------------------------------------------------
# 1. Define input files
#----------------------------------------------------------

# Get all CSV files starting with "set" in the current directory
files <- list.files(path = ".", pattern = "^set.*\\.csv$", full.names = TRUE)

# Stock genotype reference file
stock_geno <- read.csv("../stock_genotype.csv") %>% 
  select(-MitoNuc) %>% 
  mutate(Build = case_when(
    grepl("A", Stock) ~ "A",
    grepl("B", Stock) ~ "B",
    .default = "parental"
  ))

#----------------------------------------------------------
# 2. Functions
#----------------------------------------------------------

# Function to clean and parse columns from climbing data
clean_climb_data <- function(df) {
  df <- df %>%
    mutate(
      treatment = substr(Nuc, nchar(Nuc), nchar(Nuc)),
      Nuc = substr(Nuc, 1, nchar(Nuc) - 1),
      build = substr(Nuc, nchar(Nuc), nchar(Nuc)),
      Nuc = substr(Nuc, 1, nchar(Nuc) - 1),
      vial = as.integer(sapply(strsplit(vial_ID, "_"), tail, 1))
    ) %>%
    filter(vial != "all")
  return(df)
}

#----------------------------------------------------------
# 3. Read and preprocess data
#----------------------------------------------------------

# Read in all climbing data files
dfs <- lapply(files, read.csv)

# Apply column cleaning where needed
dfs <- lapply(seq_along(dfs), function(i) {
  if (all(c("Nuc", "vial_ID") %in% colnames(dfs[[i]]))) {
    clean_climb_data(dfs[[i]])
  } else {
    dfs[[i]]
  }
})

# Define set IDs (corresponds to experimental set)
sets <- c(1:7, 9, 9)

# Define columns of interest
cols <- c("Mito", "Nuc", "Sex", "rep", "slope", "build", "vial", "Set", "treatment")

# Assign set IDs and select relevant columns
dfs <- lapply(seq_along(dfs), function(i) {
  dfs[[i]] %>%
    mutate(Set = sets[i]) %>%
    select(all_of(cols))
})

#----------------------------------------------------------
# 4. Standardize naming conventions
#----------------------------------------------------------

dfs <- lapply(dfs, function(x) x %>%
                mutate(
                  Mito = str_replace(Mito, "Yak", "yak"),
                  Mito = str_replace(Mito, "Sm21", "sm21"),
                  Mito = str_replace(Mito, "Zw", "ZW")
                )
)

#----------------------------------------------------------
# 5. Combine datasets and finalize metadata
#----------------------------------------------------------

climb_df <- bind_rows(dfs) %>%
  mutate(
    treatment = case_when(
      treatment == "C" ~ "Control",
      treatment == "R" ~ "Rotenone",
      TRUE ~ treatment
    ),
    Mito = case_when(
      Mito == "Bei5" ~ "Bei05",
      TRUE ~ Mito
    )
  )

# Rename columns for consistency
colnames(climb_df) <- c("Mito", "Nuc", "Sex", "Rep", "Slope", 
                        "Build", "Vial", "Set", "Treatment")

# Merge with stock genotype reference
climb_df <- climb_df %>%
  left_join(stock_geno, by = c("Mito", "Nuc", "Build"))

#----------------------------------------------------------
# 6. Save processed data
#----------------------------------------------------------

write.csv(climb_df, "climb.csv", row.names = FALSE, quote = FALSE)