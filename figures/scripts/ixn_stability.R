library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rsample)
library(purrr)
library(broom)

weight = read.csv("weight/data/weight_adj.csv") %>%
  filter(!Mito %in% c("Ore", "375")) %>%  # remove unwanted mitochondrial lines
  mutate(Y_adj = Y_adj * 1000)%>%
  mutate(
    across(c(Mito, Treatment, Nuc, Sex, Build), as.factor)
  )

climb <- read.csv("climb/data/climb_adj.csv", header = TRUE) %>%
  filter(!Mito %in% c("375", "Ore")) %>%
  mutate(
    across(c(Mito, Treatment, Nuc, Sex, Build, Vial), as.factor)
  )

dev <- read.csv("development/data/development_adj.csv") %>%
  filter(!Mito %in% c("375", "Ore"))


## ---------------------------
## Stratify dataset by sex
## ---------------------------
climbF <- climb %>%
  filter(Sex == "F") %>%
  mutate(MitoNucBuild = paste(Mito, Nuc, Build, sep = ":"))

climbM <- climb %>%
  filter(Sex == "M") %>%
  mutate(MitoNucBuild = paste(Mito, Nuc, Build, sep = ":"))

## Stratify by sex
weightF <- weight %>% filter(Sex == "F") %>%
  na.omit()

weightM <- weight %>% filter(Sex == "M")


## ---------------------------
## Bootstrap interaction stability
## ---------------------------

# Sample n_mito Mitos without replacement, enforce a square design by keeping
# only (Treatment, Nuc, Build[, Vial]) combinations present in ALL selected Mitos,
# fit the appropriate mixed model, and return Satterthwaite ANOVA p-values.
#
# use_vial_re: include (1 | Mito:Nuc:Build:Treatment:Vial) in the model AND
#   use Vial in the combo intersection. Set FALSE when the data has a Vial
#   column that is not in the model random effects (e.g. dev data).
# covariates:  character vector of additional fixed-effect terms added after
#   Mito * Nuc * Treatment (e.g. "Larval_density").
bootstrap_mito_ixn <- function(df, n_mito = 10, n_boot = 200, seed = 42,
                                use_vial_re = "Vial" %in% names(df),
                                covariates  = NULL) {
  set.seed(seed)
  mitos <- unique(as.character(df$Mito))

  combo_cols <- if (use_vial_re) c("Treatment", "Nuc", "Build", "Vial")
                else              c("Treatment", "Nuc", "Build")

  fixed_rhs  <- paste(c("Mito * Nuc * Treatment", covariates), collapse = " + ")
  random_rhs <- if (use_vial_re)
    "(1 | Mito:Nuc:Build) + (1 | Mito:Nuc:Build:Treatment:Vial)"
  else
    "(1 | Mito:Nuc:Build)"
  model_formula <- as.formula(paste("Y_adj ~", fixed_rhs, "+", random_rhs))

  results <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    selected <- sample(mitos, n_mito, replace = FALSE)
    sub      <- df[df$Mito %in% selected, ]

    # Intersect combos across all selected Mitos to enforce a square design
    combo_sets <- lapply(selected, function(m) {
      rows <- sub[sub$Mito == m, ]
      do.call(paste, c(rows[combo_cols], list(sep = "_")))
    })
    common_combos <- Reduce(intersect, combo_sets)

    sub$combo_key <- do.call(paste, c(sub[combo_cols], list(sep = "_")))
    sub <- sub[sub$combo_key %in% common_combos, ]
    sub$combo_key <- NULL
    sub <- droplevels(sub)

    if (length(unique(sub$Mito)) < n_mito) next  # dropped a Mito; skip

    # Guard against empty Nuc × Treatment marginal cells (can arise when two
    # incomplete-data Mitos jointly cause an entire Nuc:Treatment to be removed
    # by the square intersection above — e.g. Bei54 + ZW144 in weightM both
    # lack different Build levels of Nuc=375:Rotenone, wiping it entirely).
    obs_fixed  <- nrow(unique(sub[, c("Mito", "Nuc", "Treatment")]))
    need_fixed <- length(unique(sub$Mito)) *
                  length(unique(sub$Nuc)) *
                  length(unique(sub$Treatment))
    if (obs_fixed < need_fixed) next

    fit <- tryCatch(
      suppressMessages(lmer(model_formula, data = sub, REML = TRUE)),
      error = function(e) NULL
    )

    if (is.null(fit)) next

    aov_tab <- anova(fit, ddf = "Satterthwaite")
    results[[i]] <- data.frame(
      boot    = i,
      term    = rownames(aov_tab),
      p_value = aov_tab[["Pr(>F)"]],
      F_value = aov_tab[["F value"]],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, Filter(Negate(is.null), results))
}

# Plot -log10(p) for interaction terms across bootstrap replicates.
# A dashed line marks the alpha threshold; violin + boxplot shows the distribution.
plot_ixn_stability <- function(boot_res, alpha = 0.05/4,
                               title = "Interaction stability (n = 10 Mitos/bootstrap)") {
  ixn <- boot_res[grepl(":", boot_res$term), ]

  # Proportion of bootstraps where each term is significant
  sig_prop <- ixn %>%
    group_by(term) %>%
    summarise(prop_sig = mean(p_value < alpha), .groups = "drop")

  label_df <- sig_prop %>%
    mutate(label = paste0(round(prop_sig * 100), "% sig."),
           y     = max(-log10(ixn$p_value), na.rm = TRUE) * 1.02)

  ggplot(ixn, aes(x = term, y = -log10(p_value))) +
    geom_violin(fill = "steelblue", alpha = 0.55, color = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.7, fill = "white", color = "grey30") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed",
               color = "firebrick", linewidth = 0.6) +
    geom_text(data = label_df, aes(x = term, y = y, label = label),
              size = 2.8, color = "steelblue4", vjust = 0) +
    annotate("text", x = Inf, y = -log10(alpha),
             label = paste0(" p = ", alpha), hjust = 0, vjust = -0.4,
             size = 3, color = "firebrick") +
    labs(x = NULL, y = expression(-log[10](italic(p))), title = title) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

## Run bootstraps for all five datasets
boot_res_climbF  <- bootstrap_mito_ixn(climbF,  n_mito = 10, n_boot = 1000, seed = 99)
boot_res_climbM  <- bootstrap_mito_ixn(climbM,  n_mito = 10, n_boot = 1000, seed = 99)
boot_res_weightF <- bootstrap_mito_ixn(weightF, n_mito = 10, n_boot = 1000, seed = 99)
boot_res_weightM <- bootstrap_mito_ixn(weightM, n_mito = 10, n_boot = 1000, seed = 99)
boot_res_dev     <- bootstrap_mito_ixn(dev, n_mito = 10, n_boot = 1000, seed = 99,
                                        use_vial_re = FALSE,
                                        covariates  = "Larval_density")

p1 <- plot_ixn_stability(boot_res_climbF,  title = "Female Climbing")
p2 <- plot_ixn_stability(boot_res_climbM,  title = "Male Climbing")
p3 <- plot_ixn_stability(boot_res_weightF, title = "Female Weight")
p4 <- plot_ixn_stability(boot_res_weightM, title = "Male Weight")
p5 <- plot_ixn_stability(boot_res_dev,     title = "Development Time")

final = (p1 + p2) / (p3 + p4) / (p5 + plot_spacer())
ggsave("supp_figs/ixn_bootstrap.pdf", final, width=8, height=10)
