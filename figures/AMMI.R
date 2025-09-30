############################################################
# Script: AMMI.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Perform AMMI analysis on adjusted climbing data 
#          and generate genotype-by-environment interaction 
#          biplots stratified by sex. Generates figure 2.
############################################################


## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("dplyr", "ggplot2", "agricolae", "ggrepel", "patchwork")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)

## ---------------------------------------------------------
## Load and prepare climbing dataset
## ---------------------------------------------------------
climb <- read.csv("climb/data/climb_adj.csv") %>%
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
  "parent"    = "black"
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
        grepl("sm21", name)   ~ "D.simulans",
        grepl("yak", name)    ~ "D.yakuba",
        .default              = name
      )
    )
  
  # Subset genotypes only for coloring
  mdf_gen <- mdf %>% filter(type == "GEN")
  
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
    coord_cartesian(ylim = c(-0.75, 1.1), xlim = c(-1.1, 0.9)) +
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

# Generate plots
pF <- make_plot(climbF, "Female Climb")
pM <- make_plot(climbM, "Male Climb")

# Combine plots side by side
pF + pM

# Save to PDF
ggsave("figures/climb_ammi_fig2.pdf", combined_plot, width = 10, height = 5)
