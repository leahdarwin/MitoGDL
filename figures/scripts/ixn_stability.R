###########################################################
# Script: ixn_stability.R
# Author: Leah Darwin
# Date: 2026-05-14
# Purpose: Tests stability of interaction terms via subsampling (without replacement)
#          across n_mito levels c(5, 10, 15, 20) with 1000 draws each.
############################################################

# Load libraries ----------------------------------------------------------
## -------------------------------------------------------------------
## Setup: load & install required packages
## -------------------------------------------------------------------
packages = c(
  "lme4", "lmerTest", "dplyr", "ggplot2",
  "patchwork", "ggh4x", "purrr", "broom"
)

new_pkgs = setdiff(packages, rownames(installed.packages()))
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, dependencies = TRUE)
}
invisible(lapply(packages, library, character.only = TRUE))

weight = read.csv("weight/data/weight_adj.csv") %>%
  filter(!Mito %in% c("Ore", "375")) %>%  # remove parental lines
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
  filter(!Mito %in% c("375", "Ore")) %>%
  mutate(across(c(Mito, Treatment, Nuc, Build), as.factor))


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

weightM <- weight %>% filter(Sex == "M") %>% na.omit()


## ---------------------------
## Subsampling interaction stability
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
subsample_mito_ixn <- function(df, n_mito = 10, n_draws = 200, seed = 42,
                                use_vial_re = "Vial" %in% names(df),
                                covariates  = NULL,
                                n_cores     = parallel::detectCores() - 1L) {
  mitos <- unique(as.character(df$Mito))

  combo_cols <- if (use_vial_re) c("Treatment", "Nuc", "Build", "Vial")
                else              c("Treatment", "Nuc", "Build")

  fixed_rhs  <- paste(c("Mito * Nuc * Treatment", covariates), collapse = " + ")
  random_rhs <- if (use_vial_re)
    "(1 | Mito:Nuc:Build) + (1 | Mito:Nuc:Build:Treatment:Vial)"
  else
    "(1 | Mito:Nuc:Build)"
  model_formula <- as.formula(paste("Y_adj ~", fixed_rhs, "+", random_rhs))

  one_draw <- function(i) {
    set.seed(seed + i)  # per-draw seed for reproducibility across core counts
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

    if (length(unique(sub$Mito)) < n_mito) return(NULL)

    # Guard against empty Nuc × Treatment marginal cells (can arise when two
    # incomplete-data Mitos jointly cause an entire Nuc:Treatment to be removed
    # by the square intersection above — e.g. Bei54 + ZW144 in weightM both
    # lack different Build levels of Nuc=375:Rotenone, wiping it entirely).
    obs_fixed  <- nrow(unique(sub[, c("Mito", "Nuc", "Treatment")]))
    need_fixed <- length(unique(sub$Mito)) *
                  length(unique(sub$Nuc)) *
                  length(unique(sub$Treatment))
    if (obs_fixed < need_fixed) return(NULL)

    fit <- tryCatch(
      suppressMessages(lmer(model_formula, data = sub, REML = TRUE)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    aov_tab <- anova(fit, ddf = "Satterthwaite")
    data.frame(
      draw    = i,
      term    = rownames(aov_tab),
      p_value = aov_tab[["Pr(>F)"]],
      F_value = aov_tab[["F value"]],
      stringsAsFactors = FALSE
    )
  }

  results <- parallel::mclapply(seq_len(n_draws), one_draw,
                                mc.cores = max(1L, n_cores))
  do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x) > 0, results))
}

# Plot -log10(p) for interaction terms across subsampling draws.
# Within each term, four violins (one per n_mito level) are shown side by side
# using a nested x-axis so distributions across subsampling depths are visible.
# A dashed line marks the alpha threshold.
plot_ixn_stability <- function(sub_res_list, n_mito_vals, alpha = 0.05/4,
                               title = "Interaction stability") {
  # Combine results across n_mito settings, tagging each with its n_mito value
  combined <- mapply(function(res, nm) {
    if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) return(NULL)
    res$n_mito <- nm
    res
  }, sub_res_list, n_mito_vals, SIMPLIFY = FALSE)
  combined <- do.call(rbind, Filter(Negate(is.null), combined))

  ixn <- combined[grepl(":", combined$term), ]
  ixn$n_mito_label <- factor(paste0("k = ", ixn$n_mito),
                              levels = paste0("k = ", sort(unique(ixn$n_mito))))

  # Proportion of draws where each term is significant, per n_mito level
  sig_prop <- ixn %>%
    group_by(term, n_mito_label) %>%
    summarise(prop_sig = mean(p_value < alpha), .groups = "drop") %>%
    mutate(label = as.character(round(prop_sig, 2)),
           y     = max(-log10(ixn$p_value), na.rm = TRUE) * 1.02)

  ggplot(ixn, aes(x = n_mito_label, y = -log10(p_value), fill = n_mito_label)) +
    geom_violin(alpha = 0.55, color = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.7, fill = "white", color = "grey30") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed",
               color = "firebrick", linewidth = 0.6) +
    geom_text(data = sig_prop, aes(x = n_mito_label, y = y, label = label),
              size = 2.5, color = "grey30", vjust = 0, inherit.aes = FALSE) +
    facet_wrap(~ term, scales = "free_x", strip.position = "bottom") +
    scale_fill_brewer(palette = "Blues", name = "Mitos sampled") +
    labs(x = NULL, y = expression(-log[10](italic(p))), title = title) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x      = element_text(size = 8, angle = 45, hjust = 1),
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text       = element_text(size = 9, color = "black")
    )
}

n_mito_vals <- c(5, 10, 15, 20)
n_draws     <- 1000

## Run subsamples (without replacement) for all five datasets across all n_mito settings.
## Each n_mito level runs in parallel; draws within each level are also parallelised.
## n_cores_outer controls how many n_mito levels run simultaneously; the remainder
## of available cores are split among draws within each level.
n_cores <- max(1L, parallel::detectCores() - 1L)

run_subsamples <- function(df, ...) {
  lapply(n_mito_vals, function(nm)
    subsample_mito_ixn(df, n_mito = nm, n_draws = n_draws, seed = 99,
                       n_cores = n_cores, ...))
}

sub_res_climbF  <- run_subsamples(climbF)
sub_res_climbM  <- run_subsamples(climbM)
sub_res_weightF <- run_subsamples(weightF, use_vial_re = FALSE)
sub_res_weightM <- run_subsamples(weightM, use_vial_re = FALSE)
sub_res_dev     <- run_subsamples(dev,     use_vial_re = FALSE,
                                  covariates = "Larval_density")

p1 <- plot_ixn_stability(sub_res_climbF,  n_mito_vals, title = "Female Climbing")
p2 <- plot_ixn_stability(sub_res_climbM,  n_mito_vals, title = "Male Climbing")
p3 <- plot_ixn_stability(sub_res_weightF, n_mito_vals, title = "Female Weight")
p4 <- plot_ixn_stability(sub_res_weightM, n_mito_vals, title = "Male Weight")
p5 <- plot_ixn_stability(sub_res_dev,     n_mito_vals, title = "Development Time")

final <- ((p1 + p2) / (p3 + p4) / (p5+plot_spacer())) + plot_layout(guides="collect", axis_titles="collect")
saveRDS(final, "figures/supp_figs/ixn_subsample.rds")
ggsave("figures/supp_figs/ixn_subsample.pdf", final, width = 8, height = 16)
