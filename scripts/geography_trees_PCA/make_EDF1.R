rm(list = ls())

library(dplyr)
library(ggpubr)

isotype_map<-readRDS("../../processed_data/Geo_info/isotype_map.rds")
geo_dis_plot<-readRDS("../../processed_data/genetic_similarity_and_admixutre/geo_distance_plots.RData")

combined <- ggpubr::ggarrange(ggpubr::ggarrange(isotype_map,
                                                geo_dis_plot,
                                                ncol = 1, 
                                                heights = c(0.5, 0.5)))

ggsave("../../figures/EDF1_isotype_map_geo_distance.pdf",
       plot = combined, 
       width = 7, height = 2.00995387*2, 
       units = "in", device = 'pdf')
