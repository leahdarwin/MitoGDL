############################################################
# Script: AMMI_plot.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Perform AMMI analysis on adjusted climbing data 
#          and generate genotype-by-environment interaction 
#          biplots stratified by sex. Generates figure 2.
############################################################


## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("dplyr", "ggplot2", "agricolae", "ggrepel", "patchwork", "knitr", "kableExtra")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)
source("scripts/figures/sourceData_helpers.R")

## ---------------------------------------------------------
## Load and prepare climbing dataset
## ---------------------------------------------------------
climb <- read.csv("data/climb_adj.csv") %>%
  na.omit() %>%
  group_by(Mito, Nuc, Treatment, Sex, Build) %>%
  summarise(climb = mean(Y_adj), .groups = "drop") %>%
  mutate(
    # Encode environment variables for AMMI
    NucTreat     = paste(Nuc, Treatment, sep = ":"),
    NucTreatSex  = paste(Nuc, Treatment, Sex, sep = ":")
  )

## Split by sex
climbF <- climb %>% filter(Sex == "F")
climbM <- climb %>% filter(Sex == "M")

## Define consistent color palette for mitochondrial origins
color_palette <- c(
  "Beijing"   = "#1C448E",
  "Zimbabwe"  = "#52C2BA",
  "D.yakuba"  = "#FCAB10",
  "D.simulans"= "#ED1C24",
  "parental"    = "darkgrey"
)

## ---------------------------------------------------------
## Function: Create AMMI biplot for one dataset
## ---------------------------------------------------------
make_plot <- function(df, name) {
  
  # Run AMMI model
  model <- with(df, AMMI(NucTreat, Mito, Build, climb, PC = TRUE))
  
  # Extract biplot coordinates and add mito origin labels
  mdf <- as.data.frame(model$biplot) %>%
    mutate(
      name = row.names(model$biplot),
      mitoOrig = case_when(
        grepl("B", name)      ~ "Beijing",
        grepl("Z", name)      ~ "Zimbabwe",
        grepl("siI", name)   ~ "D.simulans",
        grepl("yak", name)    ~ "D.yakuba",
        .default              = "parental"
      )
    )
  
  # Subset genotypes only for coloring
  mdf_gen <- mdf %>% filter(type == "GEN")

  # Export source data (full biplot coords: GEN + ENV rows)
  suffix <- ifelse(grepl("Female", name), "F", "M")
  save_source_data(mdf, paste0("AMMI_fig4_", suffix), "main_figs")

  # Build the biplot
  ggplot(mdf_gen, aes(x = PC1, y = PC2, color = mitoOrig)) +
    geom_point(size = 2) +
    
    # Environment vectors
    geom_segment(
      data = mdf %>% filter(type == "ENV"),
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = arrow(length = unit(0.2, "cm")),
      inherit.aes = FALSE,
      color = "black", alpha = 0.6, linewidth = 1
    ) +
    
    # Environment labels
    geom_text(
      data = mdf %>% filter(type == "ENV"),
      aes(x = PC1, y = PC2, label = name),
      inherit.aes = FALSE,
      color = "black", size = 3.2
    ) +
    
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    
    # Genotype labels (no leader lines, repel to avoid overlap)
    geom_text_repel(
      aes(label = name),
      size = 3, max.overlaps = Inf,
      box.padding = 0.2,
      segment.color = NA
    ) + 
    
    # Consistent mito colors
    scale_color_manual(
      values = setNames(color_palette[mdf_gen$mitoOrig], mdf_gen$mitoOrig),
      guide = "none"
    ) +
    
    # Styling
    theme_minimal() +
    ggtitle(name) +
    xlab(paste0("PC1 (", model$analysis$percent[1], "%)")) +
    ylab(paste0("PC2 (", model$analysis$percent[2], "%)")) +
 #   coord_cartesian(ylim = c(-0.75, 1.1), xlim = c(-1.1, 0.9)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

## ---------------------------------------------------------
## Run AMMI analysis and create plots
## ---------------------------------------------------------
mF <- with(climbF, AMMI(NucTreat, Mito, Build, climb, PC = TRUE))
print(mF$analysis)

mM <- with(climbM, AMMI(NucTreat, Mito, Build, climb, PC = TRUE))
print(mM$analysis)

## ---------------------------------------------------------
## Format p-values: fixed 4 decimal places for p >= 0.0001, switching to
## 2-significant-figure scientific notation below that (rather than the
## 0.0000 that fixed-decimal formatting would otherwise show).
## ---------------------------------------------------------
format_p <- function(p, threshold = 1e-4) {
  ifelse(p < threshold,
         formatC(signif(p, 2), format = "e", digits = 1),
         formatC(p, format = "f", digits = 4))
}

## ---------------------------------------------------------
## Gollob's F-test for IPCA axis significance (Sum Sq, Mean Sq,
## F, numerator/denominator df, p) -- printed to console as a
## LaTeX kable for copy-paste, not saved as source data.
## ---------------------------------------------------------
ipca_tab <- function(model, sex) {
  df2 <- model$ANOVA["Residuals", "Df"]  # denominator df, shared across axes
  model$analysis %>%
    mutate(
      Sex  = sex,
      IPCA = rownames(model$analysis),
      df2  = df2
    ) %>%
    select(Sex, IPCA, Df, df2, Sum.Sq, Mean.Sq, F.value, Pr.F, percent, acum)
}

ipca_sig_tab <- rbind(ipca_tab(mF, "F"), ipca_tab(mM, "M"))

print(
  kable(ipca_sig_tab %>% mutate(Pr.F = format_p(Pr.F)),
        format    = "latex",
        booktabs  = TRUE,
        linesep   = "",
        escape    = FALSE,
        digits    = c(0, 0, 0, 0, 3, 3, 2, 4, 1, 1),
        row.names = FALSE,
        col.names = c("Sex", "IPCA", "$df_1$", "$df_2$", "Sum Sq", "Mean Sq",
                      "$F$", "$p$", "\\% Var", "Cum. \\%")) %>%
    kable_styling(
      latex_options = c("hold_position"),
      full_width    = FALSE,
      font_size     = 10
    ) %>%
    column_spec(1, width = "3em")
)

# Generate plots
pF <- make_plot(climbF, "Female Climb")
pM <- make_plot(climbM, "Male Climb")

# Combine plots side by side
combined_plot <- pF + pM + plot_annotation(tag_levels = "a")

# Save to PDF
ggsave("figures/main_figs/AMMI_fig4.pdf", combined_plot, width = 8, height = 5)
