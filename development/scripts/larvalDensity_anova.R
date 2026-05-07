library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(patchwork)

dev = read.csv("development/data/development_adj.csv") 
weight = read.csv("weight/data/weight_adj.csv")%>%
  summarise(Y_adj = mean(Y_adj), .by=c(Mito,Nuc,Build,Set,Treatment,Sex))

ld_lm <- lmerTest::lmer(
  Larval_density ~ Mito * Nuc * Treatment  + (1 | Mito:Nuc:Build),
  data = dev
)
anova(ld_lm)

ld_mean = dev %>%
  summarise(Larval_density = mean(Larval_density), .by=c(Mito,Nuc,Build,Set,Treatment))

weight_F = weight %>% filter(Sex=="F") 
weight_M = weight %>% filter(Sex=="M") 

joinedF = weight_F %>%
  left_join(ld_mean, by = join_by(Mito,Nuc,Build,Set,Treatment)) %>%
  na.omit()
joinedM = weight_M %>%
  left_join(ld_mean, by = join_by(Mito,Nuc,Build,Set,Treatment)) %>%
  na.omit()

# Helper: extract R² and p-value from a simple lm as a formatted label
lm_label <- function(data) {
  fit <- lm(Y_adj * 1000 ~ Larval_density, data = data)
  s   <- summary(fit)
  r2  <- round(s$r.squared, 3)
  p   <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  paste0("R² = ", r2, "\np = ", format.pval(p, digits = 2, eps = 0.001))
}

p1 = ggplot(joinedF, aes(x=Larval_density, y=Y_adj*1000)) +
  geom_point() +
  geom_smooth(method="lm") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = lm_label(joinedF), size = 3) +
  labs(title="Females", y="Weight (mg)", x="Estimated larval density")

p2 = ggplot(joinedM, aes(x=Larval_density, y=Y_adj*1000)) +
  geom_point() +
  geom_smooth(method="lm") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = lm_label(joinedM), size = 3) +
  labs(title="Males", y="Weight (mg)", x="Estimated larval density")

p1+p2

rep_test = ld_mean %>%
  mutate(Condition = paste(Mito,Nuc,Build,sep=";")) %>%
  filter(Condition %in% c("siI;375;A","siI;Ore;A","yak;375;A","yak;Ore;A","375;375;parental","Ore;Ore;parental"))

rep_mod = lm(Larval_density ~ Condition*Treatment, data = rep_test)
anova(rep_mod)

ggplot(rep_test, aes(x =Condition, y=Larval_density, fill=Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) + # Hide outliers to avoid double-plotting
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             aes(color = Treatment)) +
  theme_linedraw() +
  scale_fill_manual(values = c("Control" = "black", "Rotenone" = "#723a83")) +
  scale_color_manual(values = c("Control" = "black", "Rotenone" = "#723a83"))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  


