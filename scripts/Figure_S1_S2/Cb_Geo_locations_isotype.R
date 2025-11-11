rm(list = ls())
getwd()

library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(maps)
library(ggpubr)
library(maps)


source("../utilities.R")

#### using the file after removal the four likely mixed isotypes
raw_data<-read.csv("../../data/20250901_Cb_2018_strains_data.csv")

indep_isotype_info<- raw_data

sample_list<- read.table("../../processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt")
nrow(sample_list)
#715

indep_isotype_info<-raw_data %>% 
  filter(strain %in% sample_list$V1) %>% 
  rename(lat=latitude) %>% 
  rename(long=longitude) %>% 
  mutate(isotype=strain)


indep_isotype_info$lat<-as.numeric(indep_isotype_info$lat)
indep_isotype_info$long<-as.numeric(indep_isotype_info$long)



#1. Hawaii
df_hw_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -176  & long < -129 & lat > 4 & lat < 46) 
hw_isotype <- as.character(df_hw_isotype$isotype)



#2. USA
df_usa_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -133  & long < -58 & lat > 24 & lat < 48) 
usa_isotype <- as.character(df_usa_isotype$isotype)

#3. Central America
df_ca_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -109.727720  & long < -77.949141 & lat > 7.149214 & lat < 24.985056) 
ca_isotype <- as.character(df_ca_isotype$isotype)

#4. South America
df_sa_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -83  & long < -33 & lat > -52 & lat < 6) 
sa_isotype <- as.character(df_sa_isotype$isotype)

#5. Europe
df_eu_isotype  <- indep_isotype_info %>%
  dplyr::filter(long > -25  & long < 46 & lat > 37 & lat < 65) 
eu_isotype <- as.character(df_eu_isotype$isotype)

#6. Africa
df_af_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -30  & long < 66 & lat > -36 & lat < 36)
af_isotype <- as.character(df_af_isotype$isotype)

#7. Australia
df_au_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 136.384136  & long < 155.751466 & lat > -39.897577 & lat < -10.910294)
au_isotype <- as.character(df_au_isotype$isotype)

#8. Taiwan
df_tw_isotype <- indep_isotype_info %>%
 dplyr::filter(long > 120  & long < 122 & lat > 21.7 & lat < 25.5)
 tw_isotype <- as.character(df_tw_isotype$isotype)

#9. Caribbean
df_car_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -67.919959  & long < -58.621015 & lat > 11.097408 & lat < 19.213168)
car_isotype <- as.character(df_car_isotype$isotype)
 


 

#10 Asia
 df_as_isotype_tmp_1 <- indep_isotype_info %>%
  dplyr::filter(long > 61.925175  & long < 155.089237 & lat > 2.977990 & lat < 53.639331) %>%
   dplyr::filter(!(long > 120  & long < 122 & lat > 21.7 & lat < 25.5))
 ### VSL2216, VSL2217, VSL2219 should be assigned as Asian samples (From India)
 ### although they don't have coordinate data
 df_as_isotype_tmp_2<-indep_isotype_info %>%
   filter(isotype %in% c("VSL2216", "VSL2217", "VSL2219"))
 df_as_isotype <- rbind( df_as_isotype_tmp_1,  df_as_isotype_tmp_2)
 as_isotype <- as.character(df_as_isotype$isotype)
 

 
 
# test - Malay archipelago
df_ma_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605)
ma_isotype <- as.character(df_ma_isotype$isotype)
length(ma_isotype)
#22samples, might be too few for creating a new group.


# test - Japan
df_jap_isotype <- indep_isotype_info %>%
  dplyr::filter(
    (long > 129.413710 & long < 147.790741 & lat > 30.091471 & lat < 45.629262) |
      (long > 123.454490 & long < 131.458193 & lat > 23.769151 & lat < 30.919667)
  )
jap_isotype <- as.character(df_jap_isotype$isotype)
length(jap_isotype)
#8samples，1 from Okinawa Islands
#might be too few to assign as a new group.


# test - New Caledonia
df_newc_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 162.481157  & long < 168.909586 & lat > -23.566620 & lat < -19.262473)
newc_isotype <- as.character(df_newc_isotype$isotype)
length(newc_isotype)
#1 sample which was previously classified as a Australian sample



# test - French Polynesia
df_fp_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -152.349878  & long < -147.559839 & lat > -18.885498 & lat < -16.161921)
fp_isotype <- as.character(df_fp_isotype$isotype)
length(fp_isotype)
#3 samples
#too few to assign as a new group.


# test - Micronesia
df_micro_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335)
micro_isotype <- as.character(df_micro_isotype$isotype)
length(micro_isotype)
#21 samples
#too few to assign as a new group.







#10 Pacific
df_pac_isotype <- indep_isotype_info %>%
  dplyr::filter(
    (long > 129.413710 & long < 147.790741 & lat > 30.091471 & lat < 45.629262) |
      (long > 123.454490 & long < 131.458193 & lat > 23.769151 & lat < 30.919667) | # Japan
      (long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605) | # Malay archipelago
      (long > 162.481157  & long < 168.909586 & lat > -23.566620 & lat < -19.262473) | # New Caledonia
      (long > -152.349878  & long < -147.559839 & lat > -18.885498 & lat < -16.161921) | # French Polynesia
      (long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335) # Micronesia
  )
pac_isotype <- as.character(df_pac_isotype$isotype)





#11. others
df_etc_isotype <- indep_isotype_info %>%
  dplyr::filter(!isotype %in% c(hw_isotype, usa_isotype, ca_isotype,car_isotype,
                               sa_isotype,eu_isotype, af_isotype,
                               au_isotype,tw_isotype, as_isotype,
                               pac_isotype))
etc_isotype <- as.character(df_etc_isotype$isotype)

## Look into the differences
etc_isotype_location <- indep_isotype_info %>% 
  dplyr::filter(isotype %in% etc_isotype)
  

  
  
  
  

## Merge                  

indep_isotype_info_geo <- indep_isotype_info %>%
  dplyr::mutate(geo = ifelse(isotype %in% hw_isotype, "Hawaii", 
                             ifelse(isotype %in% usa_isotype, "North America", 
                                    ifelse(isotype %in% ca_isotype, "Central America", 
                                           ifelse(isotype %in% sa_isotype, "South America", 
                                                  ifelse(isotype%in% eu_isotype, "Europe",
                                                         ifelse(isotype %in% af_isotype, "Africa", 
                                                                ifelse(isotype %in% au_isotype, "Australia", 
                                                                       ifelse(isotype %in% tw_isotype, "Taiwan",
                                                                              ifelse(isotype %in% car_isotype, "Caribbean",
                                                                              ifelse(isotype %in% pac_isotype, "Pacific",
                                                                                     ifelse(isotype %in% as_isotype, "Asia", "unknown"))))))))))))

indep_isotype_info_geo$lat<-as.numeric(indep_isotype_info_geo$lat)
indep_isotype_info_geo$long<-as.numeric(indep_isotype_info_geo$long)

indep_isotype_info_geo_output_tmp<-indep_isotype_info_geo %>% 
  select(isotype,lat,long,geo)


### add Cosmopolitan group
# Cosmopolitan_isotypes<-read.table("../../processed_data/Geo_info/Cb_over_100km_isotype.txt")
Cosmopolitan_isotypes<-read.table("../../processed_data/Geo_info/Cb_Cosmopolitan_isotype.txt")

indep_isotype_info_geo_output<-indep_isotype_info_geo_output_tmp %>% 
  mutate(geo = ifelse(isotype %in% Cosmopolitan_isotypes$V1,"Cosmopolitan", geo))

# write.csv(indep_isotype_info_geo_output,file="Cb_indep_isotype_info_geo.csv", row.names = FALSE, quote = FALSE)
write.csv(indep_isotype_info_geo_output,file="../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv", row.names = FALSE, quote = FALSE)

unique(indep_isotype_info_geo_output$geo)

library(RColorBrewer) #palette
display.brewer.all()


display.brewer.pal(8,"Set2")
brewer.pal(8, "Set2")
display.brewer.pal(10, "Set3")
brewer.pal(10, "Set3")
display.brewer.pal(9, "Set1")
brewer.pal(9,"Set1")






## Geographic distribution of all strains 
world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

plot_world <-ggplot2::ggplot()+ geom_map(data=world, map=world,
                                         aes(
                                           # x=long, y=lat, 
                                           map_id=region),
                                         color="gray95", fill="gray80", linewidth=0.05)+
  geom_point(data = indep_isotype_info_geo_output %>% filter(geo != "Cosmopolitan"), aes(x= long, y=lat, color = geo), shape =16, size = 1, alpha = 1.0) +
  scale_color_manual(values = geo.colours) +
  lims(x=c(-180,191),y=c(-58,84)) +
  theme_void()+
  #  theme(text = element_text(size=12)) +
  labs(color = NULL)+
  theme(legend.position = "none")+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank(),   # Conceal Tick Marks
        legend.text = element_text(size = 5))+
  theme(legend.key.size = unit(0.4, "cm"))

# ggthemes::theme_map()

plot_world


#ggsave("Cb_strains_map.pdf",plot = plot_world,width = 7.5, height = 4.37, units = "in")




### Pie charts of the strains
geo_freq <- indep_isotype_info_geo_output %>%
  dplyr::group_by(geo) %>%
  dplyr::summarize(frequency = n()) %>%
  arrange(desc(frequency))


##CHANGE THE SAVING DIRECTORIES TO RF FOR BELOW
write.csv(file = "../../processed_data/Geo_info/Cb_isotype_geo_freq.csv", geo_freq,quote = FALSE,
          row.names = FALSE)

#geo_freq$percent <- geo_freq$frequency / sum(geo_freq$frequency) * 100

plot_freq <- ggplot(geo_freq, aes(x = "", y = frequency, fill = geo)) +
 geom_bar(stat = "identity", width = 1,) +
 coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
      axis.title = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none") +
 # geom_text(aes(label = frequency), position = position_stack(vjust = 0.5)) +
 geom_text(aes (x = 1.3, label = frequency),
              color = "black",
           position = position_stack(vjust = 0.5),
           size = 3,
         hjust = 0.5) +
scale_fill_manual(values = geo.colours)+
 scale_color_manual(values = geo.colours)

plot_freq<- plot_freq +
 geom_text(aes(x = 1.65, label = geo),
           color = "black",
                      position = position_stack(vjust = 0.5),
                       #angle = 45,
                       hjust = 0.5, vjust = 0.5,
                      size = 3)
                     
plot_freq

####ALT plot freq ###### 

tmp_vec <- seq(1,nrow(geo_freq),1)
geo_freq2 <- geo_freq %>% 
  dplyr::arrange(frequency)
geo_freq2$group <-tmp_vec  
  
  

plot_freq2 <- ggplot2::ggplot(geo_freq2) +
  geom_rect(aes(xmin=0,xmax=frequency,ymin=group-0.25,ymax=group+0.25,fill=geo)) + 
  geom_text(aes(x=1,y=group+0.5,label=geo), fontface  = "bold", hjust=0, size = 1.8) +
  geom_text(aes(x=frequency+6,y=group,label=as.character(frequency)),size = 1.8) + 
  scale_fill_manual(values=geo.colours) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none')

plot_freq2





# Assemble the 2 plots
p <- ggpubr::ggarrange(ggpubr::ggarrange(plot_freq2,
                                         plot_world,
                                         ncol = 2, 
                                         labels = c("a","b"),
                                         widths = c(0.25, 0.75)))
p


saveRDS(p, file = "../../processed_data/assemble_figure_S1/isotype_map.rds")

# ggsave("Cb_isotypes_map_2025.pdf",
#        plot = p,width = 7, height = 2.00995387, units = "in",device = 'pdf',dpi=600)
# 
# ggsave("../../plots/Cb_isotypes_map_2025.pdf",
#        plot = p,width = 7, height = 2.00995387, units = "in",device = 'pdf',dpi=600)







