############################################################
# Script: weightDensity_plot.R
# Author: Leah Darwin
# Date: 2026-09-01
# Purpose: Plot weight vs. estimated larval density (Females/Males)
#          with a linear fit and R^2/p annotation. The mixed-model
#          test of Larval_density ~ Mito * Nuc * Treatment lives in
#          scripts/development/larvalDensity_anova.R; this script
#          only builds and saves the figure.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("dplyr", "ggplot2", "patchwork")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)
source("scripts/figures/sourceData_helpers.R")

## ---------------------------------------------------------
## Load and prepare data
## ---------------------------------------------------------
dev <- read.csv("data/development_adj.csv")
weight <- read.csv("data/weight_adj.csv") %>%
  summarise(Y_adj = mean(Y_adj), .by = c(Mito, Nuc, Build, Set, Treatment, Sex))

# Per-genotype mean larval density (see scripts/development/larvalDensity_anova.R
# for the mixed-model significance test on Larval_density)
ld_mean <- dev %>%
  summarise(Larval_density = mean(Larval_density), .by = c(Mito, Nuc, Build, Set, Treatment))

weight_F <- weight %>% filter(Sex == "F")
weight_M <- weight %>% filter(Sex == "M")

joined_F <- weight_F %>%
  left_join(ld_mean, by = join_by(Mito, Nuc, Build, Set, Treatment)) %>%
  na.omit()
joined_M <- weight_M %>%
  left_join(ld_mean, by = join_by(Mito, Nuc, Build, Set, Treatment)) %>%
  na.omit()

save_source_data(joined_F, "weightDensity_F", "supp_figs")
save_source_data(joined_M, "weightDensity_M", "supp_figs")

## ---------------------------------------------------------
## Helper: extract R^2 and p-value from a simple lm as a formatted label
## ---------------------------------------------------------
lm_label <- function(data) {
  fit <- lm(Y_adj * 1000 ~ Larval_density, data = data)
  s   <- summary(fit)
  r2  <- round(s$r.squared, 3)
  p   <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  paste0("R² = ", r2, "\np = ", format.pval(p, digits = 2, eps = 0.001))
}

## ---------------------------------------------------------
## Build plots (Females / Males)
## ---------------------------------------------------------
p1 <- ggplot(joined_F, aes(x = Larval_density, y = Y_adj * 1000)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = lm_label(joined_F), size = 3) +
  labs(title = "Females", y = "Weight (mg)", x = "Estimated larval density")

p2 <- ggplot(joined_M, aes(x = Larval_density, y = Y_adj * 1000)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = lm_label(joined_M), size = 3) +
  labs(title = "Males", y = "Weight (mg)", x = "Estimated larval density")

weight_ps <- p1 + p2

# Save to PDF
ggsave("figures/supp_figs/weight_density.pdf", weight_ps, width = 8, height = 3)
