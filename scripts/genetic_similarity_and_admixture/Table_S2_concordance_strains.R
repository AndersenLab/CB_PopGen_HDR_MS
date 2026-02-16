rm(list=ls())

library(dplyr)

gtcheck_raw<-read.table("../../data/concordance_raw/gtcheck.txt",
                        header = TRUE)

strain_list<-read.csv("../../data/20250901_Cb_2018_strains_data.csv") %>% 
  select(strain) %>% 
  pull()

length(strain_list)

gtcheck_TableS2<-gtcheck_raw %>% 
  filter(i %in% strain_list & j %in% strain_list)

((2018*2018)/2) - (2018/2)
## 2035153


gtcheck_TableS2<-gtcheck_TableS2 %>% 
  rename(strain1=i) %>% 
  rename(strain2=j) %>% 
  rename(`total genome-wide SNVs`=sites) %>% 
  rename(`identical alleles` = discordance)
  


write.csv(gtcheck_TableS2,
          "../../tables/TableS2.csv",
          row.names = FALSE
          )


