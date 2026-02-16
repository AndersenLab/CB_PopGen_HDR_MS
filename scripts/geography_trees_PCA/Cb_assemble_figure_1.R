rm(list = ls())

library(extrafont)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggplotify)


library(patchwork)



strain_map<-readRDS("../../processed_data/assemble_figure_1/strain_map.rds")
PCA<-readRDS("../../processed_data/assemble_figure_1/Cb_all_isotypes_PCA_plot_PC1_4.rds")
tree <- readRDS("../../processed_data/assemble_figure_1/Cb_plot_0.9_tree_equal_angle_rotated_not_rooted.rds")


# ggplot2::ggplot_build(PCA)


p1 <- strain_map + ggtitle(NULL)
p2 <- PCA       + ggtitle(NULL) + theme(legend.position = "none")
p3 <- tree      + ggtitle(NULL)


blank <- ggplot() + theme_void()

row2 <- plot_grid(
  blank, p2,
  ncol = 2,
  labels = c("c", "d"),
  rel_widths = c(1, 4)    
)





combined2 <- plot_grid(
  p1,
  row2,
  blank,
  p3,
  ncol = 1,
  label_size = 14,
  rel_heights = c(2.143951, 2.5, 3, 3.75) 
)

ggsave("../../figures/Figure1_raw_StrainMap_Tree_PCA.pdf", combined2 , width = 7, height = sum(c(2.143951, 2.5, 3, 3.75)), units = "in")





