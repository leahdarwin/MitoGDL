############################################################
# Script: climb_adj.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Apply within-set and across-set correction terms
#          to climbing assay data and save adjusted dataset.
############################################################

#----------------------------------------------------------
# 0. Package management
#----------------------------------------------------------

required_packages <- c("dplyr", "ggplot2")

# Install any missing packages
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

#----------------------------------------------------------
# 1. Load data
#----------------------------------------------------------

climb <- read.csv("climb.csv", header = TRUE)

#----------------------------------------------------------
# 2. Compute correction terms
#----------------------------------------------------------

# Within-set correction term
phenos_correct <- climb %>%
  na.omit() %>%
  group_by(Set, Sex, Treatment) %>%
  filter((Mito %in% c("yak", "sm21") & Build == "A") | Build == "parental") %>%
  summarise(Set_adj = mean(Slope), .groups = "drop")

# Across-set correction term
pos_correct <- climb %>%
  na.omit() %>%
  group_by(Treatment, Sex) %>%
  summarise(Pos_adj = mean(Slope), .groups = "drop")

#----------------------------------------------------------
# 3. Apply corrections
#----------------------------------------------------------

climb_adj <- climb %>%
  inner_join(phenos_correct, join_by(Set, Sex, Treatment)) %>%
  inner_join(pos_correct, join_by(Treatment, Sex)) %>%
  mutate(Y_adj = (Slope - Set_adj) + Pos_adj)

#----------------------------------------------------------
# 4. Save adjusted dataset
#----------------------------------------------------------

write.csv(climb_adj, "climb_adj.csv", row.names = FALSE, quote = FALSE)

