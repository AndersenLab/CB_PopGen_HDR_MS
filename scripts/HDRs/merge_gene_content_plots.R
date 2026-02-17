library(ggplot2)
library(cowplot)

p1 <- readRDS("../../processed_data/gene_diversity/HDR_I_12469000_12644000_11rg.Rds")
p2 <- readRDS("../../processed_data/gene_diversity/HDR_II_12500000_12610000_11rg.Rds")
p3 <- readRDS("../../processed_data/gene_diversity/HDR_V_838824_893260_11rg.Rds")

hdr_3reg <- cowplot::plot_grid(
  p3 + theme(legend.position = "none", axis.title.x = element_blank(),plot.title = element_text(hjust = -0.08, face = "bold", size = 14)) + ggtitle("a"),
  p2 + theme(legend.position = "none", axis.title.x = element_blank(),plot.title = element_text(hjust =-0.08, face = "bold", size = 14)) + ggtitle("b"),
  p1 + theme(legend.position = "none",plot.title = element_text(hjust = -0.08, face = "bold", size = 14)) + ggtitle("c"),
  nrow = 3,
  label_x = -0.001,
  label_y = 1.02,# move labels farther left (default is 0)
  align = "v",
  axis = "lr"
)

ggsave(plot = hdr_3reg,filename = "../../figures/EDF6_HDR_3reg_genes.png",width = 7.5,height = 8.5,device = "png",units = "in",dpi = 600)
