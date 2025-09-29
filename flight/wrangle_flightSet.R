library(dplyr)
library(stringr)

files = list.files(path=".", pattern="^flightSet.*csv$")
dfs = lapply(files,read.csv)

stock_geno = read.csv("../stock_genotype.csv") %>% 
  select(-c(MitoNuc)) %>%
  mutate(Build = case_when(grepl("A",Stock) ~ "A",
                           grepl("B",Stock) ~ "B",
                           .default = "parental"))

sets = c(1,1,2,2,3,3,4,5,6,7)
builds = c(NA, NA, "A", "B", "A", "B",NA,NA,NA,NA)
treats = c("Control", "Rotenone", rep(NA, 8))

##general edits to all dataframes
dfs = lapply(seq_along(dfs), function(i) {
  dfs[[i]] %>%
    select(c("Mito", "Nuc", "Sex", "Y")) %>%
    mutate(Set = sets[i]) %>%
    mutate(Build = builds[i]) %>%
    mutate(Treatment = treats[i]) %>% 
    filter(!grepl('test', Nuc)) %>%
    filter(!grepl('w1118', Nuc)) %>%
    filter(Y>0)
})


##changes to each individual dataframe 
for(i in c(1:2)){
  dfs[[i]] = dfs[[i]] %>%
    mutate(Build = case_when(Nuc %in% c("OreR","375") ~ "parental",
                             Nuc %in% c("OreB", "375B") ~ "B",
                             Nuc %in% c("OreA", "375A") ~ "A")) %>%
    mutate(Nuc = case_when(grepl("Ore", Nuc) ~ "Ore",
                           grepl("375", Nuc) ~ "375")) 
}

for(i in c(3:6)){
  dfs[[i]] = dfs[[i]] %>%
    mutate(Treatment = case_when(grepl("R$", Nuc) ~ "Rotenone",
                                 grepl("C$", Nuc) ~ "Control")) %>%
    mutate(Nuc = case_when(grepl("Ore", Nuc) ~ "Ore",
                           grepl("375", Nuc) ~ "375")) 
}

for(i in c(7:10)){
  dfs[[i]] = dfs[[i]] %>%
    mutate(Treatment = case_when(grepl("R$", Nuc) ~ "Rotenone",
                                 grepl("C$", Nuc) ~ "Control")) %>%
    mutate(Build = case_when(grepl("A", Nuc) ~ "A",
                             grepl("B", Nuc) ~ "B",
                             .default = "parental")) %>%
    mutate(Nuc = case_when(grepl("Ore", Nuc) ~ "Ore",
                           grepl("375", Nuc) ~ "375")) 
}


df = bind_rows(dfs)

df = df %>%
  mutate(Mito = gsub("A", "", Mito)) %>%
  mutate(Mito = gsub("B$", "", Mito)) %>%
  mutate(Mito = gsub("R$", "", Mito)) %>%
  mutate(Mito = case_when(Mito%in%c("parental","Parental") ~ Nuc,
                          .default = Mito)) %>%
  left_join(stock_geno, join_by(Mito,Nuc,Build), relationship = "many-to-one") %>%
  filter(Mito != "Bei") %>%
  mutate(Build = case_when(Mito==Nuc ~ "parental",
                           .default = Build))


write.csv(df,"flight.csv",row.names = FALSE,quote=FALSE)

