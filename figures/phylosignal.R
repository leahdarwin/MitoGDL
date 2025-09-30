############################################################
# Script: phylosignal.R
# Author: Leah Darwin
# Date: 2025-09-30
# Purpose: Test for phylogenetic signal (Blomberg’s K) in
#          fly traits across mitochondrial haplotypes.
#          Data include climbing, flight, development, and
#          weight phenotypes. 
############################################################

# Load required libraries
packages <- c("ape", "phytools", "picante", "dplyr", "tidyr", "aplot", "ggtree", "patchwork", "cowplot")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)


#-----------------------------------------------------------
# Load phylogeny 
#-----------------------------------------------------------
mttree <- read.tree("figures/GDL_mts_SR_all.raxml.bestTree")

#-----------------------------------------------------------
# Define function: calculate Blomberg’s K for a given condition
#-----------------------------------------------------------
get_k_sig <- function(condition, df) {
  # Subset dataframe to current condition
  temp_df <- df %>%
    filter(Nuc == condition["Nuc"],
           Treatment == condition["Treatment"],
           Sex == condition["Sex"])
  
  # Extract trait vector named by mitochondrial haplotype
  trait <- setNames(temp_df$Y, temp_df$Mito)
  
  #print(trait)  # Debugging: check trait vector
  
  # Run phylogenetic signal test
  invisible(
    capture.output(
      phylosig(mttree, trait, method = "K", test = TRUE)))
}



#-----------------------------------------------------------
# Load and summarize phenotypes
#-----------------------------------------------------------

## Climb
climb <- read.csv("climb/data/climb_adj.csv") %>%
  group_by(Mito, Nuc, Treatment, Sex) %>%
  na.omit() %>%
  summarise(Y = mean(Y_adj),
            se_value = sd(Y_adj) / sqrt(n()))%>%
  mutate(Mito = case_when(
    Mito == "Bei05" ~ "Bei5",
    .default = Mito
  ))

## Flight
flight <- read.csv("flight/data/flight_adj.csv") %>%
  group_by(Mito, Nuc, Treatment, Sex) %>%
  na.omit() %>%
  summarise(Y = mean(Y_adj),
            se_value = sd(Y_adj) / sqrt(n())) %>%
  mutate(Mito = case_when(
    Mito == "Bei05" ~ "Bei5",
    .default = Mito
  ))

## Development
dev <- read.csv("development/data/development_adj.csv") %>%
  group_by(Mito, Nuc, Treatment) %>%
  na.omit() %>%
  summarise(Y = mean(Y_adj),
            se_value = sd(Y_adj) / sqrt(n()),) %>%
  mutate(Sex = "MF") %>%   # No sex split for development
  mutate(Mito = case_when(
    Mito == "Bei05" ~ "Bei5",
    .default = Mito
  ))

## Weight
weight <- read.csv("weight/data/weight_adj.csv") %>%
  group_by(Mito, Nuc, Treatment, Sex) %>%
  na.omit() %>%
  summarise(Y = mean(Y_adj),
            se_value = sd(Y_adj) / sqrt(n()),) %>%
  mutate(Mito = case_when(
    Mito == "Bei05" ~ "Bei5",
    .default = Mito
  ))

#-----------------------------------------------------------
# Generate condition grids
#-----------------------------------------------------------

# For traits with sex-specific data
conditions <- expand.grid(
  Sex       = unique(climb$Sex),
  Treatment = unique(climb$Treatment),
  Nuc       = unique(climb$Nuc)
)

# For development (MF only)
cond_dev <- expand.grid(
  Sex       = unique(dev$Sex),
  Treatment = unique(climb$Treatment),
  Nuc       = unique(climb$Nuc)
)

#-----------------------------------------------------------
# Run phylogenetic signal tests
#-----------------------------------------------------------
Kclimb  <- apply(conditions, 1, get_k_sig, df = climb)
Kfly    <- apply(conditions, 1, get_k_sig, df = flight)
Kdev    <- apply(cond_dev, 1, get_k_sig, df = dev)
Kweight <- apply(conditions, 1, get_k_sig, df = weight)


#-----------------------------------------------------------
# Define a custom ggplot2 theme for tree plots
#-----------------------------------------------------------
custom_theme <- theme(
  # Backgrounds
  panel.background = element_rect(fill = "grey90", color = NA),   # Panel background
  plot.background  = element_rect(fill = "grey90", color = NA),   # Plot background
  
  # Borders and axis lines
  panel.border = element_blank(),                                 # Remove full border
  axis.line.x  = element_line(color = "black", size = 0.8),       # Keep x-axis line
  axis.line.y  = element_blank(),                                 # Remove y-axis line
  
  # Grid lines
  panel.grid.major.x = element_line(color = "black", linetype = "dotted"), # Vertical grid only
  panel.grid.major.y = element_blank(),                                   # No horizontal grid
  panel.grid.minor   = element_blank(),                                   # No minor grid
  
  # Remove y-axis details (ticks, labels, title)
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),
  
  # Legend
  legend.background = element_rect(fill = "grey90")
)

#-----------------------------------------------------------
# Tree plotting with ggtree
#-----------------------------------------------------------

# Save base plot to object for re-use
g <- ggtree(mttree) +
  geom_tiplab(align = TRUE) +
  hexpand(0.3) +
  custom_theme

# Add node labels to the tree
g + geom_text2(aes(label = node), hjust = -0.3, size = 3)

#-----------------------------------------------------------
# Create plots for each phenotype
#-----------------------------------------------------------

# Define color palette 
colors = c("Beijing" = "#1C448E","Zimbabwe" = "#52C2BA","D.simulans"="#FCAB10","D.yakuba"="#ED1C24", "375/Ore"= "darkgrey")
nuccolors = c("375"="#ef8a62","Ore"="#67a9cf")

#-----------------------------------------------------------
# Plot 1: Climbing speed (Females, Control treatment)
#-----------------------------------------------------------
p1 <- climb %>% 
  filter(Sex == "F", Treatment == "Control") %>%
  mutate(label = Mito) %>%
  ggplot(aes(x = label, y = Y, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = Y - se_value, ymax = Y + se_value),
    width = 0.2, position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = nuccolors) +
  coord_flip() +
  custom_theme +
  theme(legend.position = "none") +
  ylab("Mean climbing \nspeed (cm/s)")

#-----------------------------------------------------------
# Plot 2: Flight landing height (Females, Control treatment)
#-----------------------------------------------------------
p2 <- flight %>% 
  filter(Sex == "F", Treatment == "Control") %>%
  mutate(label = Mito) %>%
  ggplot(aes(x = label, y = Y, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = Y - se_value, ymax = Y + se_value),
    width = 0.2, position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = nuccolors) +
  coord_flip() +
  custom_theme +
  theme(legend.position = "none") +
  ylab("Mean landing height (m)")

#-----------------------------------------------------------
# Plot 3: Development time (Control treatment, MF pooled)
#-----------------------------------------------------------
p3 <- dev %>% 
  filter(Treatment == "Control") %>%
  mutate(label = Mito) %>%
  ggplot(aes(x = label, y = Y, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = Y - se_value, ymax = Y + se_value),
    width = 0.2, position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = nuccolors) +
  coord_flip() +
  custom_theme +
  theme(legend.position = "none") +
  ylab("Mean development \ntime (days)")

#-----------------------------------------------------------
# Plot 4: Weight per fly (Females, Control treatment)
#-----------------------------------------------------------
p4 <- weight %>% 
  filter(Sex == "F", Treatment == "Control") %>%
  mutate(label = Mito) %>%
  ggplot(aes(x = label, y = Y, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = Y - se_value, ymax = Y + se_value),
    width = 0.2, position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = nuccolors) +
  coord_flip() +
  custom_theme +
  theme(legend.position = "bottom") +
  ylab("Mean weight (mg)")

#-----------------------------------------------------------
# Combine plots with tree and spacers
#-----------------------------------------------------------
combined_plot <- p1 %>%
  insert_left(g, width = 1.2) %>%
  insert_right(plot_spacer(), width = 0.05) %>%
  insert_right(p2) %>%
  insert_right(plot_spacer(), width = 0.05) %>%
  insert_right(p3) %>%
  insert_right(plot_spacer(), width = 0.05) %>%
  insert_right(p4)

# Save to PDF
ggsave("figures/tree_pheno_fig4.pdf", combined_plot, width = 14, height = 5)
