library(aplot)
library(ggplot2)
library(ggtree)
library(cowplot)
library(patchwork)

setwd("/Users/leahdarwin/Documents/drand/gdlPhenotypeExperiment")

tree = read.tree("GDL_mts_SR_all.raxml.bestTree")

climb = read.csv("dist/climb/climb_adj.csv")
flight = read.csv("dist/flight/flight.csv")
dev = read.csv("dist/development/development_adj.csv")
weight = read.csv("dist/weight/weight_adj.csv")

colors = c("Beijing" = "#1C448E","Zimbabwe" = "#52C2BA","D.simulans"="#FCAB10","D.yakuba"="#ED1C24", "375/Ore"= "darkgrey")
nuccolors = c("375"="#ef8a62","Ore"="#67a9cf")

custom_theme <- theme(
  panel.background = element_rect(fill = "grey90", color = NA),  # Grey background
  plot.background = element_rect(fill = "grey90", color = NA),   # Grey outer background
  panel.border = element_blank(),  # Fully remove full border
  axis.line.x = element_line(color = "black", size = 0.8),  # Keep dark y-axis border
  axis.line.y = element_blank(),  # Remove x-axis line
  panel.grid.major.x = element_line(color = "black", linetype = "dotted"),  # Only vertical grid lines
  panel.grid.major.y = element_blank(),  # No horizontal grid lines
  panel.grid.minor = element_blank(),  # No minor grid lines
  axis.text.y = element_blank(),  # Remove y-axis tick labels
  axis.ticks.y = element_blank(),  # Remove y-axis ticks
  axis.title.y = element_blank(),  # Remove y-axis title
  legend.background = element_rect(fill = "grey90")
)

ggtree(tree, branch.length = "none")+ geom_tiplab(align=TRUE)+ hexpand(.3)

g = ggtree(tree)+ geom_tiplab(align=TRUE) + hexpand(.3)
g +geom_text2(aes(label = node), hjust = -0.3, size = 3)
g %>%collapse(node=36) %>%collapse(node=29)


p1 = climb %>% 
  filter(Sex == "F", treatment=="C") %>%
  group_by(Mito,Nuc) %>%
  summarise(
    mean_value = mean(slope),
    se_value = sd(slope) / sqrt(n()),  # Standard Error
    .groups = "drop"
  ) %>%
  mutate(label=Mito) %>%
  ggplot(aes(x = label, y = mean_value, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +  # Bar plot
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2, position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values=nuccolors)+
  coord_flip() + 
  custom_theme +
  theme(legend.position='none') +
  ylab("Mean climbing \nspeed (cm/s)")
p1

p2 = flight %>% 
  filter(Sex == "F", Treatment=="Control") %>%
  mutate(Mito = case_when(Mito=="Bei05"~"Bei5",
                          .default = Mito)) %>%
  group_by(Mito,Nuc) %>%
  summarise(
    mean_value = mean(Y),
    se_value = sd(Y) / sqrt(n()),  # Standard Error
    .groups = "drop"
  ) %>%
  mutate(label=Mito) %>%
  ggplot(aes(x = label, y = mean_value, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +  # Bar plot
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2, position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values=nuccolors)+
  coord_flip() + 
  custom_theme +
  theme(legend.position='none') +
  ylab("Mean landing height (m)") 
p2


p3 = dev %>% 
  filter( Treatment=="Control") %>%
  group_by(Mito,Nuc) %>%
  mutate(Mito = case_when(Mito=="Bei05"~"Bei5",
                          Mito=="zim53"~"Zim53",
                          .default = Mito)) %>%
  summarise(
    mean_value = mean(TSE),
    se_value = sd(TSE) / sqrt(n()),  # Standard Error
    .groups = "drop"
  ) %>%
  mutate(label=Mito) %>%
  ggplot(aes(x = label, y = mean_value, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +  # Bar plot
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2, position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values=nuccolors)+
  coord_flip() + 
  custom_theme +
  theme(legend.position='none') +
  ylab("Mean development \ntime (days)") 
p3


p4 = weight %>% 
  filter(Sex == "F", Treatment=="Control") %>%
  mutate(Mito = case_when(Mito=="Bei05"~"Bei5",
                          Mito=="zim53"~"Zim53",
                          .default = Mito)) %>%
  group_by(Mito,Nuc) %>%
  summarise(
    mean_value = mean(weight_per_fly*1000),
    se_value = sd(weight_per_fly*1000) / sqrt(n()),  # Standard Error
    .groups = "drop"
  ) %>%
  mutate(label=Mito) %>%
  ggplot(aes(x = label, y = mean_value, fill = Nuc)) +
  geom_col(position = position_dodge(width = 0.8)) +  # Bar plot
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                width = 0.2, position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values=nuccolors)+
  coord_flip() + 
  custom_theme +
  theme(legend.position='bottom') +
  ylab("Mean weight (mg)") 
p4 

p1  %>% insert_left(g, width=1.2) %>% insert_right(plot_spacer(), width = 0.05)%>%insert_right(p2) %>%insert_right(plot_spacer(), width = 0.05) %>% insert_right(p3) %>%insert_right(plot_spacer(), width = 0.05) %>% insert_right(p4)
