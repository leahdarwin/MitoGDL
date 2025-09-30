############################################################
# Script: mothersCurse_plot.R
# Author: Leah Darwin
# Date: 2025-09-29
# Purpose: Analyze fly trait datasets (climb, flight, weight)
#          to assess variance, coefficient of variation (CV),
#          sex differences, and correlations. Generates
#          bootstrapped plots, variance/correlation tables,
#          and combined summary figures.
############################################################

# Load libraries ----------------------------------------------------------
packages <- c("dplyr", "ggplot2", "kableExtra", 
              "ggpubr", "patchwork", "ggtext", "tidyr")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

set.seed(105)


# Helper functions --------------------------------------------------------

# Compute coefficient of variation
cv <- function(x) sd(x) / mean(x)

# Bootstrap CV values
boot_cv <- function(x, n = 500) {
  replicate(n, cv(sample(x, replace = TRUE)))
}

# Extract variance test results (M vs F)
get_var_df <- function(df) {
  df %>%
    group_by(Nuc, Treatment) %>%
    summarise(
      var_test = list(var.test(
        Y_adj[Sex == "M"],
        Y_adj[Sex == "F"],
        alternative = "greater"
      )),
      .groups = "drop"
    ) %>%
    mutate(
      F_stat = sapply(var_test, \(x) x$statistic),
      p_value = sapply(var_test, \(x) x$p.value)
    ) %>%
    select(Nuc, Treatment, F_stat, p_value)
}

# Make bootstrapped CV plot
make_cv_plot <- function(df, trait, trait_color) {
  boot_climb <- df %>% 
    group_by(Mito, Nuc, Sex, Treatment) %>%
    summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
    group_by(Sex, Nuc, Treatment) %>%
    summarise(boot = list(boot_cv(Y_adj)), .groups = "drop") %>%
    tidyr::unnest(boot)
  
  ggplot(boot_climb, aes(
    y = boot,
    x = interaction(Nuc, Treatment, sep = "\n"),
    fill = Sex, color = Sex
  )) +
    geom_point(position = position_jitterdodge(jitter.width = 0.3),
               size = 0.6, alpha = 0.5) +
    geom_boxplot(color = "black", alpha = 0.7, outliers = FALSE) +
    labs(title = trait, x = "", y = "Bootstrapped CV") +
    theme_bw() +
    scale_color_manual(values = sex_colors) +
    scale_fill_manual(values = sex_colors) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(size = 0.5, colour = "black"),
          axis.line.y = element_line(size = 0.5, colour = "black")) +
    theme(
      plot.title = element_textbox_simple(
        size = 16, lineheight = 1, halign = 0.5,
        fill = trait_color, padding = margin(5, 5, 5, 5),
        margin = margin(0, 0, 5, 47)
      ),
      plot.title.position = "plot",
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.9),
      text = element_text(size = 16)
    )
}

# Bootstrap-based p-values
get_boot_pvals <- function(df, n = 1000) {
  boot_climb <- df %>%
    group_by(Mito, Nuc, Sex, Treatment) %>%
    summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
    group_by(Sex, Nuc, Treatment) %>%
    summarise(boot = list(boot_cv(Y_adj, n)), .groups = "drop") %>%
    tidyr::unnest(boot)
  
  sexbootclimb <- boot_climb %>%
    group_by(Nuc, Treatment) %>%
    tidyr::pivot_wider(names_from = Sex, values_from = boot, values_fn = list)
  
  sexbootclimbgeq <- Map(function(a, b) Map(`>=`, a, b), sexbootclimb$`F`, sexbootclimb$M)
  
  means <- boot_climb %>%
    group_by(Sex, Nuc, Treatment) %>%
    summarise(boot = mean(boot), .groups = "drop")
  
  tibble(
    Nuc = sexbootclimb$Nuc,
    Treatment = sexbootclimb$Treatment,
    fMean = means %>% filter(Sex == "F") %>% pull(boot),
    mMean = means %>% filter(Sex == "M") %>% pull(boot),
    pval = sapply(sexbootclimbgeq, function(x) sum(unlist(x)) / n)
  )
}

# Correlation plot and test
make_corr_plot <- function(df) {
  mvf <- df %>%
    filter(Sex %in% c("F", "M")) %>%
    tidyr::pivot_wider(names_from = Sex, values_from = Y_adj)
  
  ggplot(mvf, aes(x = F, y = M)) +
    geom_point(size = 3)
}

get_corr_df = function(df){ 
  
  maleF = df %>% filter(Sex=="F") 
  maleT = df %>% filter(Sex=="M") 
  
  mvf = maleF %>% 
    left_join(maleT, join_by(Mito,Nuc,Treatment)) 
  
  mvf = mvf %>% 
    group_by(Nuc, Treatment) %>% 
    summarise( cor_test = list(cor.test(Y_adj.x, Y_adj.y, method = "pearson")), .groups = "drop" ) %>% 
    mutate( r = sapply(cor_test, function(x) x$estimate), p = sapply(cor_test, function(x) x$p.value) ) %>% 
    select(-cor_test) 
}

# Define sex color palette
sex_colors <- c("F" = "lightgrey", "M" = "gray34")

# Load datasets -----------------------------------------------------------

climb <- read.csv("climb/data/climb_adj.csv") %>%
  na.omit() %>%
  group_by(Sex, Treatment, Nuc, Mito, Build) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
  group_by(Sex, Treatment, Nuc, Mito) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
  group_by(Sex, Treatment, Nuc) %>%
  mutate(group_mean = mean(Y_adj),
         group_var = var(Y_adj),
         group_cv = group_var / group_mean)

flight <- read.csv("flight/data/flight_adj.csv") %>%
  na.omit() %>%
  group_by(Sex, Treatment, Nuc, Mito, Build) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
  group_by(Sex, Treatment, Nuc, Mito) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
  group_by(Sex, Treatment, Nuc) %>%
  mutate(group_mean = mean(Y_adj),
         group_var = var(Y_adj),
         group_cv = group_var / group_mean)

weight <- read.csv("weight/data/weight_adj.csv") %>%
  na.omit() %>%
  mutate(Y_adj = Y_adj * 1000) %>%
  group_by(Sex, Treatment, Nuc, Mito, Build) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
  group_by(Sex, Treatment, Nuc, Mito) %>%
  summarise(Y_adj = mean(Y_adj), .groups = "drop_last") %>%
  group_by(Sex, Treatment, Nuc) %>%
  mutate(group_mean = mean(Y_adj),
         group_var = var(Y_adj)) %>%
  filter(!(Mito == "Bei54" & Treatment == "Rotenone"),
         !(Mito == "ZW144" & Treatment == "Rotenone" & Nuc == "Ore")) %>%
  mutate(group_cv = group_var / group_mean)

# Run analyses ------------------------------------------------------------

# Bootstrapped CV plots + p-values
pbootclimb   <- make_cv_plot(climb,  "Climbing Velocity", "#b3cde3")
pbootweight  <- make_cv_plot(weight, "Weight", "#fbb4ae")
pbootflight  <- make_cv_plot(flight, "Flight Performance", "#ccebc5")

bootpvals_climb  <- get_boot_pvals(climb)  %>% mutate(Trait = "climb")
bootpvals_weight <- get_boot_pvals(weight) %>% mutate(Trait = "weight")
bootpvals_flight <- get_boot_pvals(flight) %>% mutate(Trait = "flight")

cvboot_df <- rbind(bootpvals_climb, bootpvals_weight, bootpvals_flight)

# Combine bootstrapped plots
combined_plot <- pbootclimb + theme(legend.position = "none") +
  pbootweight + theme(legend.position = "none") + labs(y = "") +
  pbootflight + labs(y = "")

# Save to PDF
ggsave("figures/mothers_curse_fig3.pdf", combined_plot, width = 15, height = 5)

# Correlation results
corrdf <- rbind(get_corr_df(weight) %>% mutate(Trait="weight"),
                get_corr_df(flight) %>% mutate(Trait="flight"),
                get_corr_df(climb)  %>% mutate(Trait="climb")) %>%
  as.data.frame()

# Variance test results
vardf <- rbind(get_var_df(weight) %>% mutate(Trait="weight"),
               get_var_df(flight) %>% mutate(Trait="flight"),
               get_var_df(climb)  %>% mutate(Trait="climb"))

# Merge correlation + variance + CV
corvar <- corrdf %>%
  inner_join(vardf,   join_by(Nuc, Treatment, Trait)) %>%
  inner_join(rbind(
    climb  %>% select(Nuc, Treatment, Sex, group_var, group_cv) %>% mutate(Trait = "climb"),
    weight %>% select(Nuc, Treatment, Sex, group_var, group_cv) %>% mutate(Trait = "weight"),
    flight %>% select(Nuc, Treatment, Sex, group_var, group_cv) %>% mutate(Trait = "flight")
  ), join_by(Nuc, Treatment, Trait), relationship = "many-to-many") %>%
  mutate(r = paste0(round(r,4), " (", round(p,4), ")"),
         F_stat = paste0(round(F_stat,4), " (", round(p_value,4), ")")) %>%
  select(-c(p, p_value)) %>%
  select(Trait, Nuc, Treatment, Sex, group_var, F_stat, r)

# Table output
kable(cvboot_df, format = "latex", digits = 4) %>%
  kable_styling()
