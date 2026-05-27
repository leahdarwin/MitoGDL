############################################################
# Script: phenoCorr_plot.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Correlation analysis of phenotypic traits
#          (weight, climb, flight, development) for
#          males and females across nuclear backgrounds.
############################################################

# Load libraries ----------------------------------------------------------
packages <- c("dplyr", "ggplot2", "kableExtra", 
              "ggpubr", "patchwork", "ggtext", 
              "tidyr", "ggfortify", "ggsci", "stringr", "knitr")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)


# -------------------------------------------------------------------
# Load datasets
# -------------------------------------------------------------------
weight <- read.csv("data/weight.csv")
climb  <- read.csv("data/climb.csv")
dev    <- read.csv("data/development.csv")
flight <- read.csv("data/flight.csv")

# -------------------------------------------------------------------
# Helper function: summarize trait by group
# -------------------------------------------------------------------
summarize_trait <- function(df, trait_col, new_name, scale_factor = 1, sex_filter = NULL) {
  df %>%
    na.omit() %>%
    { if (!is.null(sex_filter)) filter(., Sex == sex_filter) else . } %>%
    group_by(Mito, Nuc, Treatment) %>%
    summarise(!!new_name := mean(.data[[trait_col]] * scale_factor),
              .groups = "drop")
}

# -------------------------------------------------------------------
# Summarize each trait by sex
# -------------------------------------------------------------------
weight_F <- summarize_trait(weight, "weight_per_fly", "weight_F", scale_factor = 1000, sex_filter = "F")
weight_M <- summarize_trait(weight, "weight_per_fly", "weight_M", scale_factor = 1000, sex_filter = "M")

climb_F  <- summarize_trait(climb, "Slope", "climb_F", sex_filter = "F")
climb_M  <- summarize_trait(climb, "Slope", "climb_M", sex_filter = "M")

flight_F <- summarize_trait(flight, "Y", "flight_F", sex_filter = "F")
flight_M <- summarize_trait(flight, "Y", "flight_M", sex_filter = "M")

dev_MF   <- summarize_trait(dev, "TSE", "dev")

# -------------------------------------------------------------------
# Merge trait summaries (by sex)
# -------------------------------------------------------------------
merged_F <- dev_MF %>%
  full_join(climb_F,  by = c("Mito","Nuc","Treatment")) %>%
  full_join(flight_F, by = c("Mito","Nuc","Treatment")) %>%
  full_join(weight_F, by = c("Mito","Nuc","Treatment")) %>%
  na.omit()

merged_M <- dev_MF %>%
  full_join(weight_M, by = c("Mito","Nuc","Treatment")) %>%
  full_join(climb_M,  by = c("Mito","Nuc","Treatment")) %>%
  full_join(flight_M, by = c("Mito","Nuc","Treatment")) %>%
  na.omit()

# -------------------------------------------------------------------
# Clean up mito names for plotting
# -------------------------------------------------------------------
standardize_mito <- function(df) {
  df %>%
    mutate(Mito_orig = case_when(
      grepl("Bei", Mito) ~ "Bei",
      grepl("Z",   Mito) ~ "Zim",
      TRUE               ~ Mito
    ))
}

merged_F <- standardize_mito(merged_F)
merged_M <- standardize_mito(merged_M)

# -------------------------------------------------------------------
# Phenotype combinations
# -------------------------------------------------------------------
combs_F <- combn(colnames(merged_F)[4:7], 2, simplify = FALSE)
combs_M <- lapply(combs_F, function(x) gsub("_F", "_M", x))

# -------------------------------------------------------------------
# Axis label mapper
# -------------------------------------------------------------------
get_axisLab <- function(pheno) {
  case_when(
    grepl("climb",  pheno) ~ "Climbing velocity (cm/s)",
    grepl("flight", pheno) ~ "Mean landing height (m)",
    grepl("dev",    pheno) ~ "Development time (days)",
    grepl("weight", pheno) ~ "Weight (mg)"
  )
}

# -------------------------------------------------------------------
# Plotting function
# -------------------------------------------------------------------
make_plot <- function(comb, data) {
  ggplot(data, aes_string(x = comb[1], y = comb[2], color = "Nuc")) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Ore" = "#67a9cf", "375" = "#ef8a62")) +
    geom_smooth(data = subset(data, Nuc == "375"), method = lm, formula = y ~ x) +
    geom_smooth(data = subset(data, Nuc == "Ore"), method = lm, formula = y ~ x) +
    xlab(get_axisLab(comb[1])) +
    ylab(get_axisLab(comb[2])) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_bw() +
    labs(color = "Nuclear Genotype") +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# -------------------------------------------------------------------
# Correlation testing
# -------------------------------------------------------------------
get_cor <- function(comb, nuc, data) {
  col1 <- data %>% filter(Nuc == nuc) %>% pull(comb[1])
  col2 <- data %>% filter(Nuc == nuc) %>% pull(comb[2])
  test <- cor.test(col1, col2)
  return(c(comb[1], comb[2], nuc, test$estimate, test$conf.int[1], test$conf.int[2], test$p.value))
}

# -------------------------------------------------------------------
# Run correlation tests and print tables; return plots
# -------------------------------------------------------------------
run_corr_tests <- function(combs, merged_df) {
  corr_results <- rbind(
    lapply(combs, get_cor, "375", merged_df),
    lapply(combs, get_cor, "Ore", merged_df)
  )

  corr_tab <- as.data.frame(do.call(rbind, corr_results)) %>%
    `colnames<-`(c("Phenotype 1", "Phenotype 2", "Nuc", "Correlation coeff", "CI lower", "CI upper", "p-value")) %>%
    mutate(
      `Phenotype 1` = get_axisLab(`Phenotype 1`),
      `Phenotype 2` = get_axisLab(`Phenotype 2`),
      `Correlation coeff` = as.numeric(`Correlation coeff`),
      `CI lower`          = as.numeric(`CI lower`),
      `CI upper`          = as.numeric(`CI upper`),
      `p-value`           = as.numeric(`p-value`)
    )

  # Print nicely formatted LaTeX table
  print(
    kable(corr_tab,
          format    = "latex",
          booktabs  = TRUE,
          linesep   = "",
          escape    = FALSE,
          digits    = rep(3, 7),
          col.names = c("Phenotype 1", "Phenotype 2", "Nuc", "$r$", "$r_{lower}$", "$r_{upper}$", "$p$")) %>%
      kable_styling(
        latex_options = c("hold_position"),
        full_width    = FALSE,
        font_size     = 10
      ) %>%
      column_spec(2, width = "10em")
  )

  # Return plots for external assembly
  lapply(combs, make_plot, merged_df)
}

# -------------------------------------------------------------------
# Run for Females and Males
# -------------------------------------------------------------------
plots_F <- run_corr_tests(combs_F, merged_F)
plots_M <- run_corr_tests(combs_M, merged_M)

# -------------------------------------------------------------------
# Assemble into a single figure using patchwork
# -------------------------------------------------------------------
female_panel <- wrap_plots(plots_F, nrow = 2, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Females") &
  theme(legend.position = "bottom")

male_panel <- wrap_plots(plots_M, nrow = 2, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Males") &
  theme(legend.position = "bottom")

combined <- wrap_elements(female_panel) / wrap_elements(male_panel) +
  plot_annotation(tag_levels = "a")

ggsave(
  "figures/supp_figs/phenoCorr_figS2.pdf",
  combined,
  width = 8, height = 10
)

