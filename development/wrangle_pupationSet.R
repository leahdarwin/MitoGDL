library("car")
library("ggplot2")
library("dplyr")
library("stringr")

stock_geno = read.csv("../stock_genotype.csv")

files = list.files(path=".", pattern="^pupationSet.*csv$")
dfs = lapply(files,read.csv)

edit_cols = function(df){
  df$treatment = substr(df$Nuc, nchar(df$Nuc), nchar(df$Nuc))
  df$Nuc = substr(df$Nuc, 1, nchar(df$Nuc) - 1)
  df$build = substr(df$Nuc, nchar(df$Nuc), nchar(df$Nuc))
  df$Nuc = substr(df$Nuc, 1, nchar(df$Nuc) - 1)
  df$vial = sapply(strsplit(df$vial_ID, "_"), tail, 1)
  df = subset(df,vial!="all")
}

sets = c(2:9,9)

dfs = lapply(seq_along(dfs), function(i) {
  dfs[[i]] %>%
    mutate(Set = sets[i])
})

day_weights = c(6.5:12.5)

wrangle = function(df){
  
  data = df %>%
    rowwise() %>%
    mutate(TSE = sum(c_across(3:9)*day_weights)/sum(c_across(3:9))) %>%
    mutate(build = case_when(
      str_detect(Stock, "A") ~ "A",
      str_detect(Stock, "B") ~ "B",
      TRUE ~ "parental"
    )) %>% 
    mutate(total = sum(c_across(3:9))) %>%
    select(-c(3:9))
  
  data = inner_join(data,stock_geno, by="Stock")
  
}

wrangled_dfs = lapply(dfs, wrangle)
all_sets = bind_rows(wrangled_dfs)

colnames(all_sets) = c("Stock", "Vial", "Treatment", "Set", "TSE", "Build", "Total", "MitoNuc", "Mito", "Nuc")

all_sets = all_sets %>%
  mutate(Mito = case_when(Mito == "zim53" ~ "Zim53",
                          .default = Mito))


write.csv(all_sets, file="development.csv", row.names = FALSE, quote=FALSE)


