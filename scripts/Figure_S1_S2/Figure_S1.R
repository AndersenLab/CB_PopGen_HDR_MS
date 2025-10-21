rm(list = ls())

library(dplyr)

isotype_map<-readRDS("../../processed_data/assemble_figure_S1/isotype_map.rds")
isotype_map

geo_dis_plot<-readRDS("../../processed_data/assemble_figure_S1/geo_distance_plots.RData")
geo_dis_plot


library(ggpubr)

combined <- ggpubr::ggarrange(ggpubr::ggarrange(isotype_map,
                                                geo_dis_plot,
                                                ncol = 1, 
                                                # labels = c("c","d"),
                                                heights = c(0.5, 0.5)))
combined

ggsave("Figure_S1.pdf",
       plot = combined, 
       width = 7, height = 2.00995387*2, 
       units = "in", device = 'pdf')



