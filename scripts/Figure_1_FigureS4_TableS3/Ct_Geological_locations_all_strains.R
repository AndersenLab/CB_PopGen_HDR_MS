rm(list = ls())
getwd()

list.files()

library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(maps)
library(ggpubr)
library(maps)


source("../utilities.R")



# raw_data <- read.table("../../data/Ct_spreadsheet.txt",
#                    sep = '\t',header = TRUE)
raw_data<-read.csv("../../data/20250627_c_tropicalis_strain_data.csv")


indep_strain_info<- raw_data


indep_strain_info<-raw_data %>% 
  rename(lat=latitude) %>% 
  rename(long=longitude) %>% 
  select(strain,isotype,lat,long)


indep_strain_info$lat<-as.numeric(indep_strain_info$lat)
indep_strain_info$long<-as.numeric(indep_strain_info$long)
 


#1. Central America
df_ca_strain <- indep_strain_info %>%
  dplyr::filter(long > -97.8  & long < -76 & lat > 9 & lat < 10.5) 
ca_strain <- as.character(df_ca_strain$isotype)

#2. South America
df_sa_strain <- indep_strain_info %>%
  dplyr::filter(long > -72  & long < -52 & lat > -13 & lat < 6) 
sa_strain <- as.character(df_sa_strain$isotype)

#3. Hawaii
df_hw_strain <- indep_strain_info %>%
  dplyr::filter(long > -176  & long < -129 & lat > 4 & lat < 46)  
hw_strain <- as.character(df_hw_strain$isotype)

#4. New zealand
df_nz_strain <- indep_strain_info %>%
  dplyr::filter(long > 165  & long < 180) %>%
  dplyr::filter(lat > -50 & lat < -25)
nz_strain <- as.character(df_nz_strain$isotype)


#5. Australia
df_au_strain <- indep_strain_info %>%
  dplyr::filter(long > 110  & long < 160) %>%
  dplyr::filter(lat > -50 & lat < -10) 
au_strain <- as.character(df_au_strain$isotype)

#6. Caribbean
df_car_strain <- indep_strain_info %>%
  dplyr::filter(long > -80  & long < -50) %>%
  dplyr::filter(lat > 10 & lat < 20) 
car_strain <- as.character(df_car_strain$isotype)


#7. Africa
df_af_strain <- indep_strain_info %>%
  dplyr::filter(long > -30  & long < 70) %>%
  dplyr::filter(lat > -40 & lat < 30) 
af_strain <- as.character(df_af_strain$isotype)

#8. Europe
df_eu_strain <- indep_strain_info %>%
  dplyr::filter(long > -10  & long < 25) %>%
  dplyr::filter(lat > 35 & lat < 57) 
eu_strain <- as.character(df_eu_strain$isotype)


#9. Taiwan
df_tw_strain <- indep_strain_info %>%
  dplyr::filter(long > 50  & long < 150) %>%
  dplyr::filter(lat > 12 & lat < 75) 
tw_strain <- as.character(df_tw_strain$isotype)


# #10. Pacific
# df_pac_strain <- indep_strain_info %>%
#   dplyr::filter(
#     (long > 129.413710 & long < 147.790741 & lat > 30.091471 & lat < 45.629262) |
#       (long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605) | # Malay Archipelago
#       (long > 162.481157  & long < 168.909586 & lat > -23.566620 & lat < -19.262473) | # New Caledonia
#       (long > -152.349878  & long < -147.559839 & lat > -18.885498 & lat < -16.161921) | # French Polynesia
#       (long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335) # Micronesia
#   )
# pac_strain <- as.character(df_pac_strain$isotype)
# 




#11. Malay Archipelago
df_mal_strain <- indep_strain_info %>%
  dplyr::filter(long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605)  # Malay Archipelago
mal_strain <- as.character(df_mal_strain$isotype)



#12. Micronesia
df_mic_strain <- indep_strain_info %>%
  dplyr::filter(long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335) # Micronesia
mic_strain <- as.character(df_mic_strain$isotype)








#0. others
df_etc_strain <- indep_strain_info %>%
  dplyr::filter(!isotype %in% c(af_strain, au_strain, car_strain,
                                eu_strain,hw_strain, ca_strain, 
                                nz_strain,sa_strain, tw_strain,
                                #### pac_strain,
                                mal_strain, mic_strain
  ))
etc_strain <- as.character(df_etc_strain$isotype)



### Merge                  
indep_strain_info_geo <- indep_strain_info %>%
  dplyr::mutate(geo = ifelse(isotype %in% hw_strain, "Hawaii", 
                             ifelse(isotype %in% ca_strain, "Central America", 
                                    ifelse(isotype %in% au_strain, "Australia", 
                                           ifelse(isotype %in% af_strain, "Africa", 
                                                  ifelse(isotype %in% car_strain, "Caribbean",
                                                         ifelse(isotype %in% eu_strain, "Europe", 
                                                                ifelse(isotype %in% sa_strain, "South America", 
                                                                       ifelse(isotype %in% nz_strain, "New Zealand", 
                                                                              ifelse(isotype %in% tw_strain, "Taiwan", 
                                                                                     # ifelse(isotype %in% pac_strain, "Pacific", 
                                                                                     ifelse(isotype %in% mal_strain,"Malay Archipelago",
                                                                                            ifelse(isotype %in% mic_strain,"Micronesia","Unknown"))))))))))))
unique(indep_strain_info_geo$geo)


######## Look into the differences #######
etc_strain_location <- indep_strain_info %>% 
  dplyr::filter(strain %in% etc_strain)
etc_strain_location
# write.csv(etc_strain_location,"undefined_region_strains.csv", quote = FALSE,row.names = FALSE)
  
  


library(RColorBrewer) #palette
display.brewer.all()


display.brewer.pal(8,"Set2")
brewer.pal(8, "Set2")
display.brewer.pal(10, "Set3")
brewer.pal(10, "Set3")
display.brewer.pal(9, "Set1")
brewer.pal(9,"Set1")



# geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
#                  "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
#                  "Taiwan" = "#E5C494", "Unknown" = "black","Europe" = "red")



# 
# geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
#                  "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F",
#                  "Taiwan" = "#E5C494", "North America" = "#A65628", "Europe" = "#377EB8",
#                  "Asia" = "#E41A1C", "Pacific" = "#984EA3","unknown" = 'grey')

## Geographic distribution of all strains 
world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

plot_world <-ggplot2::ggplot()+ geom_map(data=world, map=world,
                                         aes(x=long, y=lat, map_id=region),
                                         color="gray95", fill="gray80", linewidth=0.05)+
  geom_point(data = indep_strain_info_geo, aes(x= long, y=lat, color = geo), shape =16, size = 1.2, alpha = 1.0) +
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
geo_freq <- indep_strain_info_geo %>%
  dplyr::group_by(geo) %>%
  dplyr::summarize(frequency = n()) %>%
  arrange(desc(frequency))


write.csv(file = "Ct_strains_geo_freq.csv", geo_freq,quote = FALSE,
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
  
  
  
  # dplyr::mutate(tmp=gsub("\n",".",geo)) %>%
  # dplyr::mutate(tmp2=gsub("Central.","C. ",tmp)) %>%
  # dplyr::mutate(geo2=gsub("South.","S. ",tmp)) %>%
  # dplyr::select(-tmp,-tmp2,-geo) %>%
  # dplyr::rename(geo=geo2)

plot_freq2 <- ggplot2::ggplot(geo_freq2) +
  geom_rect(aes(xmin=0,xmax=frequency,ymin=group-0.25,ymax=group+0.25,fill=geo)) + 
  geom_text(aes(x=1,y=group+0.42,label=geo),
            fontface = "bold", hjust=0, size = 1.8) +
  geom_text(aes(x=frequency+15,y=group,label=as.character(frequency)),size = 1.8) + 
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

#save image 4*4 inches




# Assemble the 2 plots
p <- ggpubr::ggarrange(ggpubr::ggarrange(plot_freq2,
                                         plot_world,
                                         ncol = 2, 
                                         labels = c("a","b"),
                                         widths = c(0.25, 0.75)))
p


ggsave("Ct_strains_map_2025.pdf",
       plot = p,width = 7.5, height = 2.153522, units = "in",device = 'pdf',dpi=600)

saveRDS(p, file = "../../processed_data/assemble_figure_1/strain_map.rds")


# The ratio of width to height for input data (the world map without Antarctica) is 2.612. 
# If the width is 7.5 inches, the wide ratio of map and frequency pot is 3:1 (75% map, 25% frequency), 
# the height of the map should be 7.5*0.75 / 2.612 = 2.153522. 







