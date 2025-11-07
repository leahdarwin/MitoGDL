############################################################
# Script: buildCorr_plot.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Calculate correlations between builds for each trait.
############################################################

# -------------------------------------------------------------------
# Load libraries
# -------------------------------------------------------------------
packages <- c("dplyr", "ggplot2", "patchwork")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

# -------------------------------------------------------------------
# Load and preprocess input data for each trait
# -------------------------------------------------------------------

# Each CSV file contains adjusted phenotype data (Y_adj)
climb  <- read.csv("climb/data/climb_adj.csv")

# Add a Sex column to development data (if not present)
dev    <- read.csv("development/data/development_adj.csv") %>% 
  mutate(Sex = "MF")

weight <- read.csv("weight/data/weight_adj.csv") %>% mutate(Y_adj = Y_adj*1000)
flight <- read.csv("flight/data/flight_adj.csv") 


# -------------------------------------------------------------------
# Define plotting function
# -------------------------------------------------------------------

make_plot <- function(df, pheno) {
  
  # Summarize phenotype means per genotype and condition
  df_p <- df %>%
    na.omit() %>%
    group_by(Mito, Nuc, Sex, Build, Treatment) %>%
    summarise(Y_adj = mean(Y_adj)) %>%
    filter(Build != "parental")  # Exclude parental lines
  
  # Merge Build A and Build B data for paired comparison
  build <- df_p %>% 
    filter(Build == "A") %>%
    full_join(df_p %>% filter(Build == "B"),
              by = join_by(Mito, Nuc, Sex, Treatment)) %>%
    na.omit()  # Remove incomplete pairs (only A or only B)
  
  # Compute correlation between Build A and Build B
  ct <- cor.test(build$Y_adj.x, build$Y_adj.y)
  
  # Extract R² and p-value
  r2 <- ct$estimate^2
  pval <- ct$p.value
  
  # Create scatterplot with annotation and 1:1 reference line
  p <- ggplot(build, aes(x = Y_adj.x, y = Y_adj.y)) +
    geom_point() +
    annotate("text", 
             x = Inf, y = -Inf, 
             label = paste0("R² = ", round(r2, 3),
                            "\np = ", signif(pval, 3)),
             hjust = 1.1, vjust = -0.5,
             size = 5,
             color="darkblue") +
    geom_abline(intercept = 0, slope = 1, 
                color = "red", linetype = "dotted") +
    theme_bw() +
    labs(x = "Build A", y = "Build B", title = pheno)
  
  return(p)
}


# -------------------------------------------------------------------
# Generate plots for all traits
# -------------------------------------------------------------------

# Bundle data frames and labels into a list
traits = list(
  list(climb,  "Climbing"),
  list(dev,    "Development"),
  list(flight, "Flight"),
  list(weight, "Weight")
)

# Apply plotting function to each trait
plots = lapply(traits, function(x) make_plot(x[[1]], x[[2]]))

# Combine all plots into a single layout using patchwork
final = wrap_plots(plots)

# Save to PDF
ggsave("figures/supp_figs/buildCorr_figS1.pdf", final, width = 7, height = 6)
