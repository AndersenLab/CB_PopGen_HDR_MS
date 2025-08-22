rm(list = ls())
library(dplyr)


isotype_by_lineage_from_Nic_raw<-read.table("../../data/From_Nic/isotype_byLineage_GeoLocAdmCol.tsv",
         header = TRUE,
         sep = '\t',
         comment.char  = "")


Cb_indep_isotype_info_geo<-read.csv("../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv")




### generate updated dataset
isotype_by_lineage<-isotype_by_lineage_from_Nic_raw %>% 
  select(isotype,Sublineage,Lineage) %>% 
  left_join(Cb_indep_isotype_info_geo,by = "isotype") 
# %>% 
#   mutate(Lineage = ifelse(Lineage == "TH",
#                           "TH",
#                           Lineage))



unique(isotype_by_lineage$Lineage)
unique(isotype_by_lineage$Sublineage)





#########################################
######### Temperate lineage #############
#########################################
Temperate_lineage<-isotype_by_lineage %>% 
  filter(Lineage == "Temperate")


write.csv(Temperate_lineage,
          "../../processed_data/Temperate_hard_filtered_and_LD_pruned/Temperate_lineage.csv",
          row.names = FALSE,
          quote = FALSE
)



#########################################
######### Tropical lineage #############
#########################################
Tropical_lineage<-isotype_by_lineage %>% 
  filter(Lineage == "Tropical")


write.csv(Tropical_lineage,
          "../../processed_data/Tropical_hard_filtered_and_LD_pruned/Tropical_lineage.csv",
          row.names = FALSE,
          quote = FALSE
)












#############################################
######## Taiwan isotype_by_lineage ##########
#############################################
Taiwan_samples<-isotype_by_lineage %>% 
  filter(geo == "Taiwan") %>% 
  group_by(Lineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are five lineages in the Taiwan geo groups
### TWD2 and TWD3 lineages have too few isotypes.
### Hence, calculate TD1, Tropical, and TH


TD1_tw<-isotype_by_lineage %>% 
  filter(geo == "Taiwan") %>% 
  filter(Lineage == "TD1")

Tropical_tw<-isotype_by_lineage %>% 
  filter(geo == "Taiwan") %>% 
  filter(Lineage == "Tropical") 

TH_tw<-isotype_by_lineage %>% 
  filter(geo == "Taiwan") %>% 
  filter(Lineage == "TH")


Taiwan_lineage<-rbind(TD1_tw,Tropical_tw,TH_tw)


write.csv(Taiwan_lineage,
          "../../processed_data/pi_theta_d_by_lineage/Taiwan_lineage/Taiwan_lineage.csv",
          row.names = FALSE,
          quote = FALSE
          )







#############################################
######## Australia isotype_by_lineage ######
#############################################
Australia_Lineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Australia")%>% 
  group_by(Lineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are two lineages in the Australia geo groups
### AD and Tropical
  
Australia_Sublineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Australia")%>% 
  group_by(Sublineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are three sublineages in the Australia geo groups
### AD, Tropical_Islandic and Tropical_Tropical



AD_au<-isotype_by_lineage %>% 
  filter(geo == "Australia") %>% 
  filter(Lineage == "AD")

Tropical_au_lineage<-isotype_by_lineage %>% 
  filter(geo == "Australia") %>% 
  filter(Lineage == "Tropical") 

Tropical_Islandic_au_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Australia") %>% 
  filter(Sublineage == "Tropical_Islandic") 

Tropical_Tropical_au_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Australia") %>% 
  filter(Sublineage == "Tropical_Tropical") 


Australia_lineage<-rbind(AD_au,Tropical_au_lineage)
Australia_sublineage<-rbind(Tropical_Islandic_au_sublineage,Tropical_Tropical_au_sublineage)




write.csv(Australia_lineage,
          "../../processed_data/pi_theta_d_by_lineage/Australia_lineage/Australia_lineage.csv",
          row.names = FALSE,
          quote = FALSE
)

write.csv(Australia_sublineage,
          "../../processed_data/pi_theta_d_by_lineage/Australia_sublineage/Australia_sublineage.csv",
          row.names = FALSE,
          quote = FALSE
)















#############################################
######## Pacific isotype_by_lineage ######
#############################################
Pacific_Lineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Pacific")%>% 
  group_by(Lineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are three lineages in the Pacific geo groups
### Indonesia, Temperate, and Tropical
### Indonesia and Temperate have too few sample.
### Hence, calculate Tropical only


Pacific_Sublineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Pacific")%>% 
  group_by(Sublineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are five sublineages in the Pacific geo groups of Tropical Lineage
### Indonesia, Temperate,Temperate_Adjacent, Temperate_Adjacent, and Tropical_Tropical. 
### Indonesia, Temperate, Temperate_Adjacent, and Tropical_Islandic have too few sample.
### Hence, calculate Tropical_Tropical only



Tropical_pa<-isotype_by_lineage %>% 
  filter(geo == "Pacific") %>% 
  filter(Lineage == "Tropical")

Tropical_Tropical_pa_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Pacific") %>% 
  filter(Sublineage == "Tropical_Tropical") 

Pacific_lineage<-rbind(Tropical_pa)
Pacific_sublineage<-rbind(Tropical_Tropical_pa_sublineage)




write.csv(Pacific_lineage,
          "../../processed_data/pi_theta_d_by_lineage/Pacific_lineage/Pacific_lineage.csv",
          row.names = FALSE,
          quote = FALSE
)

write.csv(Pacific_sublineage,
          "../../processed_data/pi_theta_d_by_lineage/Pacific_sublineage/Pacific_sublineage.csv",
          row.names = FALSE,
          quote = FALSE
)














#############################################
######## Asia isotype_by_lineage ######
#############################################
Asia_Lineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Asia")%>% 
  group_by(Lineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are Four lineages in the Asia geo groups
### Hubei, Kerala, TH, and Tropical
### Hubei and TH have too few sample.
### Hence, calculate Kerala and Tropical only


Asia_Sublineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Asia")%>% 
  group_by(Sublineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are three sublineages in the Asia geo groups of Tropical Lineage
### Temperate_Adjacent, Tropical_Islandic, and Tropical_Tropical. 
### Tropical_Islandic have too few sample.
### Hence, calculate Temperate_Adjacent and Tropical_Tropical only



Kerala_as<-isotype_by_lineage %>% 
  filter(geo == "Asia") %>% 
  filter(Lineage == "Kerala")
Tropical_as<-isotype_by_lineage %>% 
  filter(geo == "Asia") %>% 
  filter(Lineage == "Tropical")


Temperate_Adjacent_as_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Asia") %>% 
  filter(Sublineage == "Temperate_Adjacent")
Tropical_Tropical_as_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Asia") %>% 
  filter(Sublineage == "Tropical_Tropical") 

Asia_lineage<-rbind(Kerala_as,Tropical_as)
Asia_sublineage<-rbind(Temperate_Adjacent_as_sublineage,
                       Tropical_Tropical_as_sublineage)




write.csv(Asia_lineage,
          "../../processed_data/pi_theta_d_by_lineage/Asia_lineage/Asia_lineage.csv",
          row.names = FALSE,
          quote = FALSE
)

write.csv(Asia_sublineage,
          "../../processed_data/pi_theta_d_by_lineage/Asia_sublineage/Asia_sublineage.csv",
          row.names = FALSE,
          quote = FALSE
)













#############################################
######## Central America isotype_by_lineage ######
#############################################
Central_America_Lineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Central America")%>% 
  group_by(Lineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are two lineages in the Central America geo groups
### TH, and Tropical
### TH, have only one isotype.
### For now, don't re-calculate it


Central_America_Sublineage_samples<-isotype_by_lineage %>% 
  filter(geo == "Central America")%>% 
  group_by(Sublineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are two sublineages in the Central America geo groups of Tropical Lineage
### Tropical_Central, and Tropical_Tropical. 
### Tropical_Central have 100 isotypes, Tropical_Tropical have 20 isotypes.
### Hence, calculate both 

Tropical_Central_ca_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Central America") %>% 
  filter(Sublineage == "Tropical_Central")
Tropical_Tropical_ca_sublineage<-isotype_by_lineage %>% 
  filter(geo == "Central America") %>% 
  filter(Sublineage == "Tropical_Tropical") 


Central_America_sublineage<-rbind(Tropical_Central_ca_sublineage,
                                  Tropical_Tropical_ca_sublineage)


write.csv(Central_America_sublineage,
          "../../processed_data/pi_theta_d_by_lineage/Central_America_sublineage/Central_America_sublineage.csv",
          row.names = FALSE,
          quote = FALSE
)













#############################################
######## South America isotype_by_lineage ######
#############################################
South_America_Lineage_samples<-isotype_by_lineage %>% 
  filter(geo == "South America")%>% 
  group_by(Lineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are two lineages in the South America geo groups
### TH, and Tropical
### TH, have only one isotype.
### For now, don't re-calculate it


South_America_Sublineage_samples<-isotype_by_lineage %>% 
  filter(geo == "South America")%>% 
  group_by(Sublineage) %>% 
  summarise(Frequency = n()) %>% 
  ungroup()
### There are three sublineages in the South America geo groups of Tropical Lineage
### Temperate_Adjacent, Tropical_Central, and Tropical_Tropical. 
### Temperate_Adjacenonly have 8 isotypes
### Hence, calculate both Tropical_Central, and Tropical_Tropical. 



Tropical_Central_sa_sublineage<-isotype_by_lineage %>% 
  filter(geo == "South America") %>% 
  filter(Sublineage == "Tropical_Central")
Tropical_Tropical_sa_sublineage<-isotype_by_lineage %>% 
  filter(geo == "South America") %>% 
  filter(Sublineage == "Tropical_Tropical") 


South_America_sublineage<-rbind(Tropical_Central_sa_sublineage,
                                  Tropical_Tropical_sa_sublineage)



write.csv(South_America_sublineage,
          "../../processed_data/pi_theta_d_by_lineage/South_America_sublineage/South_America_sublineage.csv",
          row.names = FALSE,
          quote = FALSE
)


