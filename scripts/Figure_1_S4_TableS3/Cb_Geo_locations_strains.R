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

# raw_data<-read_excel("../../data/C. briggsae WI strain info_20250507.xlsx")
# raw_data<-read.csv("../../data/20250626_c_briggsae_strain_data.csv")
raw_data<-read.csv("../../data/20250901_Cb_2018_strains_data.csv")

indep_strain_info<- raw_data



indep_strain_info$lat<-as.numeric(indep_strain_info$lat)
indep_strain_info$long<-as.numeric(indep_strain_info$long)





#1. Hawaii
df_hw_strain <- indep_strain_info %>%
  dplyr::filter(long > -176  & long < -129 & lat > 4 & lat < 46) 
hw_strain <- as.character(df_hw_strain$strain)



#2. USA
df_usa_strain <- indep_strain_info %>%
  dplyr::filter(long > -133  & long < -58 & lat > 24 & lat < 48) 
usa_strain <- as.character(df_usa_strain$strain)

#3. Central America
df_ca_strain <- indep_strain_info %>%
  dplyr::filter(long > -109.727720  & long < -77.949141 & lat > 7.149214 & lat < 24.985056) 
ca_strain <- as.character(df_ca_strain$strain)

#4. South America
df_sa_strain <- indep_strain_info %>%
  dplyr::filter(long > -83  & long < -33 & lat > -52 & lat < 6) 
sa_strain <- as.character(df_sa_strain$strain)

#5. Europe
df_eu_strain  <- indep_strain_info %>%
  dplyr::filter(long > -25  & long < 46 & lat > 37 & lat < 65) 
eu_strain <- as.character(df_eu_strain$strain)

#6. Africa
df_af_strain <- indep_strain_info %>%
  dplyr::filter(long > -30  & long < 66 & lat > -36 & lat < 36)
af_strain <- as.character(df_af_strain$strain)

#7. Australia
df_au_strain <- indep_strain_info %>%
  dplyr::filter(long > 136.384136  & long < 155.751466 & lat > -39.897577 & lat < -10.910294)
au_strain <- as.character(df_au_strain$strain)

#8. Taiwan
df_tw_strain <- indep_strain_info %>%
 dplyr::filter(long > 120  & long < 122 & lat > 21.7 & lat < 25.5)
 tw_strain <- as.character(df_tw_strain$strain)

#9. Caribbean
df_car_strain <- indep_strain_info %>%
  dplyr::filter(long > -67.919959  & long < -58.621015 & lat > 11.097408 & lat < 19.213168)
car_strain <- as.character(df_car_strain$strain)
 


#10 Asia
 df_as_strain_tmp_1 <- indep_strain_info %>%
  dplyr::filter(long > 61.925175  & long < 155.089237 & lat > 2.977990 & lat < 53.639331) %>%
   dplyr::filter(!(long > 120  & long < 122 & lat > 21.7 & lat < 25.5))
 ### VSL2216, VSL2217, VSL2219 should be assigned as Asian samples (From India)
 ### although they don't have coordinate data
 df_as_strain_tmp_2<-indep_strain_info %>%
   filter(strain %in% c("VSL2216", "VSL2217", "VSL2219"))
 df_as_strain <- rbind( df_as_strain_tmp_1,  df_as_strain_tmp_2)
 as_strain <- as.character(df_as_strain$strain)
 

 
 
# test - Malay archipelago
df_ma_strain <- indep_strain_info %>%
  dplyr::filter(long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605)
ma_strain <- as.character(df_ma_strain$strain)
length(ma_strain)
#31

# test - Japan
df_jap_strain <- indep_strain_info %>%
  dplyr::filter(
    (long > 129.413710 & long < 147.790741 & lat > 30.091471 & lat < 45.629262) |
      (long > 123.454490 & long < 131.458193 & lat > 23.769151 & lat < 30.919667)
  )
jap_strain <- as.character(df_jap_strain$strain)
length(jap_strain)
#13


# test - New Caledonia
df_newc_strain <- indep_strain_info %>%
  dplyr::filter(long > 162.481157  & long < 168.909586 & lat > -23.566620 & lat < -19.262473)
newc_strain <- as.character(df_newc_strain$strain)
length(newc_strain)
#3 

# test - French Polynesia
df_fp_strain <- indep_strain_info %>%
  dplyr::filter(long > -152.349878  & long < -147.559839 & lat > -18.885498 & lat < -16.161921)
fp_strain <- as.character(df_fp_strain$strain)
length(fp_strain)
#3 

# test - Micronesia
df_micro_strain <- indep_strain_info %>%
  dplyr::filter(long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335)
micro_strain <- as.character(df_micro_strain$strain)
length(micro_strain)
#86 samples







#10 Pacific
df_pac_strain <- indep_strain_info %>%
  dplyr::filter(
    (long > 129.413710 & long < 147.790741 & lat > 30.091471 & lat < 45.629262) |
      (long > 123.454490 & long < 131.458193 & lat > 23.769151 & lat < 30.919667) | # Japan
      (long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605) | # Malay archipelago
      (long > 162.481157  & long < 168.909586 & lat > -23.566620 & lat < -19.262473) | # New Caledonia
      (long > -152.349878  & long < -147.559839 & lat > -18.885498 & lat < -16.161921) | # French Polynesia
      (long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335) # Micronesia
  )
pac_strain <- as.character(df_pac_strain$strain)





#11. others
df_etc_strain <- indep_strain_info %>%
  dplyr::filter(!strain %in% c(hw_strain, usa_strain, ca_strain,car_strain,
                               sa_strain,eu_strain, af_strain,
                               au_strain,tw_strain, as_strain,
                               pac_strain))
etc_strain <- as.character(df_etc_strain$strain)

## Look into the differences
etc_strain_location <- indep_strain_info %>% 
  dplyr::filter(strain %in% etc_strain)
  

  
  
  
  

## Merge                  

indep_strain_info_geo <- indep_strain_info %>%
  dplyr::mutate(geo = ifelse(strain %in% hw_strain, "Hawaii", 
                             ifelse(strain %in% usa_strain, "North America", 
                                    ifelse(strain %in% ca_strain, "Central America", 
                                           ifelse(strain %in% sa_strain, "South America", 
                                                  ifelse(strain%in% eu_strain, "Europe",
                                                         ifelse(strain %in% af_strain, "Africa", 
                                                                ifelse(strain %in% au_strain, "Australia", 
                                                                       ifelse(strain %in% tw_strain, "Taiwan",
                                                                              ifelse(strain %in% car_strain, "Caribbean",
                                                                              ifelse(strain %in% pac_strain, "Pacific",
                                                                                     ifelse(strain %in% as_strain, "Asia", "unknown"))))))))))))

indep_strain_info_geo$lat<-as.numeric(indep_strain_info_geo$lat)
indep_strain_info_geo$long<-as.numeric(indep_strain_info_geo$long)

indep_strain_info_geo_output_tmp<-indep_strain_info_geo %>% 
  select(strain,isotype,lat,long,geo)


### export strains geo info for each strain ##
indep_info_geo_for_each_strain<-indep_strain_info_geo_output_tmp %>% 
  rename(strain_geo=geo)
write.csv(indep_info_geo_for_each_strain,file="indep_info_geo_for_each_strain.csv", row.names = FALSE, quote = FALSE)
write.csv(indep_info_geo_for_each_strain,file="../../processed_data/Geo_info/indep_info_geo_for_each_strain.csv", row.names = FALSE, quote = FALSE)




### add Cosmopolitan group
# Cosmopolitan_isotypes<-read.table("../../processed_data/Geo_info/Cb_over_100km_isotype.txt")
Cosmopolitan_isotypes<-read.table("../../processed_data/Geo_info/Cb_Cosmopolitan_isotype.txt")
indep_strain_info_geo_output<-indep_strain_info_geo_output_tmp %>% 
  mutate(geo = ifelse(isotype %in% Cosmopolitan_isotypes$V1,"Cosmopolitan", geo))



write.csv(indep_strain_info_geo_output,file="Cb_indep_strain_info_geo.csv", row.names = FALSE, quote = FALSE)
write.csv(indep_strain_info_geo_output,file="../../processed_data/Geo_info/Cb_indep_strain_info_geo.csv", row.names = FALSE, quote = FALSE)

unique(indep_strain_info_geo_output$geo)

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
  geom_point(data = indep_info_geo_for_each_strain, aes(x= long, y=lat, color = strain_geo), shape =16, size = 1, alpha = 1.0) +
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




### Pie charts of the strains
geo_freq <- indep_info_geo_for_each_strain %>%
  dplyr::group_by(strain_geo) %>%
  dplyr::summarize(frequency = n()) %>%
  arrange(desc(frequency))


##CHANGE THE SAVING DIRECTORIES TO RF FOR BELOW
write.csv(file = "strain_geo_freq.csv", geo_freq,quote = FALSE,
          row.names = FALSE)

#geo_freq$percent <- geo_freq$frequency / sum(geo_freq$frequency) * 100

plot_freq <- ggplot(geo_freq, aes(x = "", y = frequency, fill = strain_geo)) +
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
 geom_text(aes(x = 1.65, label = strain_geo),
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
  geom_rect(aes(xmin=0,xmax=frequency,ymin=group-0.2,ymax=group+0.2,fill=strain_geo)) + 
  geom_text(aes(x=1,y=group+0.48,label=strain_geo),
            fontface = "bold", hjust=0, size = 2.5) +
  geom_text(aes(x=frequency+40,y=group,label=as.character(frequency)),size = 2.5) + 
  scale_fill_manual(values=geo.colours) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none')+
  xlim(0, max(geo_freq2$frequency) + 50)


plot_freq2





# Assemble the 2 plots
p <- ggpubr::ggarrange(ggpubr::ggarrange(plot_freq2,
                                         plot_world,
                                         ncol = 2, 
                                         labels = c("a","b"),
                                         widths = c(0.2, 0.8)))
p


# ggsave("Cb_strains_map_2025.pdf",
#        plot = p,width = 7.5, height = 2.153522, units = "in",device = 'pdf',dpi=600)


# The ratio of width to height for input data (the world map without Antarctica) is 2.612.
# 7/5*4 * 2.612
ggsave("Cb_strains_map_2025.pdf",
       plot = p,width = 7, height = 2.143951, units = "in",device = 'pdf',dpi=600)


saveRDS(p, file = "../../processed_data/assemble_figure_1/strain_map.rds")




