
rm(list=ls())

library(readr)
library(dplyr)

lineage<- read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byLineage_GeoLocAdmCol_20250909.tsv") %>% 
  select(isotype,Lineage)


Tropical_AD<-lineage %>% 
  filter(Lineage %in% c("Tropical","AD"))

Tropical_KD<-lineage %>% 
  filter(Lineage %in% c("Tropical","KD"))

Tropical_TD1<-lineage %>% 
  filter(Lineage %in% c("Tropical","TD1"))

Tropical_Temperate<-lineage %>% 
  filter(Lineage %in% c("Tropical","Temperate"))

Tropical_TH<-lineage %>% 
  filter(Lineage %in% c("Tropical","TH"))



# write.table(Tropical_AD,"../../processed_data/Dxy_Tropical/Tropical_AD.txt",
#             quote = FALSE,  
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = '\t')
# 
# 
# 
# write.table(Tropical_KD,"../../processed_data/Dxy_Tropical/Tropical_KD.txt",
#             quote = FALSE,  
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = '\t')
# 
# 
# 
# write.table(Tropical_TD1,"../../processed_data/Dxy_Tropical/Tropical_TD1.txt",
#             quote = FALSE,  
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = '\t')
# 
# 
# write.table(Tropical_Temperate,"../../processed_data/Dxy_Tropical/Tropical_Temperate.txt",
#             quote = FALSE,  
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = '\t')
# 
# 
# 
# 
# write.table(Tropical_TH,"../../processed_data/Dxy_Tropical/Tropical_TH.txt",
#             quote = FALSE,  
#             row.names = FALSE,
#             col.names = FALSE,
#             sep = '\t')


diversity_and_divergence



write.table(Tropical_AD,"../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_AD.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')



write.table(Tropical_KD,"../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_KD.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')



write.table(Tropical_TD1,"../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_TD1.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')


write.table(Tropical_Temperate,"../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_Temperate.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')




write.table(Tropical_TH,"../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_TH.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')



