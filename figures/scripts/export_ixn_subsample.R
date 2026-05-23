library(ggplot2)
library(patchwork)

final <- readRDS("../supp_figs/ixn_subsample.rds") 

ggsave("../supp_figs/ixn_subsample.pdf", final, width = 10, height = 9.5)
