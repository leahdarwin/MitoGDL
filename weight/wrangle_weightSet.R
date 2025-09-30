library(dplyr)

stock_geno = read.csv("../stock_genotype.csv")

files = list.files(path=".", pattern="^weightSet.*csv$")
dfs = lapply(files,read.csv)

weight_df = bind_rows(dfs)
weight_df = merge(weight_df, stock_geno, all=TRUE)
weight_df = na.omit(weight_df)

weight_df$weight_per_fly = (weight_df$Weight_flies-weight_df$Weight)/weight_df$Inds
weight_df = subset(weight_df,weight_per_fly>0)
weight_df$Treatment = gsub("Rotenone ", "Rotenone", weight_df$Treatment)

weight_df$Build <- ifelse(grepl("A", weight_df$Stock), "A",
                          ifelse(grepl("B", weight_df$Stock), "B",
                                 "parental"))

weight_df = weight_df %>%
  mutate(Mito = case_when(Mito == "zim53" ~ "Zim53",
                          .default = Mito))

write.csv(weight_df, file="weight.csv", row.names = FALSE, quote=FALSE)
