
rm(list=ls())

library(ggtree)
library(ggplot2)
library(treeio)
library(ape)
library(dplyr)
library(ggnewscale)
library(phytools)

library(RColorBrewer)

source("../utilities.R")



display.brewer.all()
display.brewer.pal(11,"RdYlBu")
brewer.pal(11,"RdYlBu")
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")

isotype_map_0.6_raw<-treeio::read.tree("../../processed_data/LD_pruned_trees/LD_0.6/phy_file_LD_0.6.phy.contree")
isotype_map_0.7_raw<-treeio::read.tree("../../processed_data/LD_pruned_trees/LD_0.7/phy_file_LD_0.7.phy.contree")
isotype_map_0.8_raw<-treeio::read.tree("../../processed_data/LD_pruned_trees/LD_0.8/phy_file_LD_0.8.phy.contree")
isotype_map_0.9_raw<-treeio::read.tree("../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy.contree")

### set midpoint rooting of phylogenetic tree
isotype_map_0.6<-phytools::midpoint_root(isotype_map_0.6_raw)
isotype_map_0.7<-phytools::midpoint_root(isotype_map_0.7_raw)
isotype_map_0.8<-phytools::midpoint_root(isotype_map_0.8_raw)
isotype_map_0.9<-phytools::midpoint_root(isotype_map_0.9_raw)

## write rooted tree
write.tree(isotype_map_0.9,
           file = "../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy.contree.rooted")



##### annotation_maps
annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)






########### new function ###########
############ plot the tree with  heatmaps #########

# function: plot the tree
plot_tree_heatmaps <- function(tree_file, 
                                    annotation_map_file,
                                    annotation_temperature_file,
                                    # annotation_Thomas_file,
                                    # annotation_MAF_file,
                                    tree_layout="circular",
                                    add_tiplab=TRUE,
                                    heatmap_lat_offset = 0.0013,
                                    # heatmap_temp_offset = 0.0007,
                                    # heatmap_geo_offset = 0.00055,
                                    heatmap_geo_offset = 0.0024
                                    # heatmap_thomas_offset = 0.0040,
                                    # heatmap_MAF_offset = 0.0051
                               ){
  
  data_tree_file<-ggtree::fortify(tree_file)
  annotation_maps<- annotation_map_file
  # annotation_temperature <- annotation_temperature_file
  # annotation_Thomas<-annotation_Thomas_file
  # annotation_MAF<-annotation_MAF_file
  
  heatmap<- annotation_maps%>%
    select(lat)
  row.names(heatmap)<-annotation_maps$isotype
  heatmap$lat<-abs(heatmap$lat)
  
  # heatmap_2<- annotation_temperature%>%
  #   select(average_temp)
  # row.names(heatmap_2)<-annotation_temperature$isotype
  
  heatmap_3<- annotation_maps%>%
    select(geo)
  row.names(heatmap_3)<-annotation_maps$isotype
  
  
  # plot tree
  p0<-ggtree::ggtree(data_tree_file, aes(col=geo), layout=tree_layout, size=0.2) %<+% annotation_maps +
    geom_tippoint(aes(color=geo), size=0.1, alpha =0.6, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
    
    # geom_tiplab(aes(label=label, col=geo), size = 0.5, hjust=-0.5, align=TRUE, linesize=0.1, alpha=0.6)+
    # geom_tiplab(aes(label=NA, col=NA), size=0.8)+
    theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  if (add_tiplab==TRUE) {
    p0<- p0+geom_tiplab(aes(label=label),
                        color = "grey50",
                        linetype = "dashed",
                        size = 0.6, hjust = -0.5,
                        align = TRUE, linesize = 0.05,
                        alpha = 0.6)
  }
  
  
  # p1<- ggtree::open_tree(p0, 90) %>% ggtree::rotate_tree(90) 
  p1<- ggtree::open_tree(p0, 1) %>% ggtree::rotate_tree(90)
  
  
  # lat_color<-brewer.pal(11,"RdYlBu")[2:10]
  p2<-p1+new_scale_fill()+new_scale_color()
  p3<-gheatmap(p2, heatmap, offset = heatmap_lat_offset, width=0.05, font.size=2.5, 
               color = NULL, hjust = -0.1, colnames_level=colnames(heatmap), 
               colnames_angle=0, 
               legend_title=" Absolute latitude",
               custom_column_labels  = rep("", ncol(heatmap))
  ) + 
    scale_fill_gradient2(
      # low = "#D73027", mid = "#FFFFBF",high = "#4575B4",
      low = "#D73027", mid = "#FFFFBF",high = "#4575B4",
      midpoint = mean(c(max(heatmap$lat, na.rm = TRUE), min(heatmap$lat, na.rm = TRUE))),
      # midpoint = mean(c(max(heatmap$lat), min(heatmap$lat))),
      na.value = "white",  # don't show NA in the heatmap
      name = "Absolute latitude"
    )
  
  # p4<-p3+new_scale_fill()+new_scale_color()
 
  
  p6<-p3 + new_scale_fill() + new_scale_color()
  p7<-gheatmap(p6, heatmap_3, offset = heatmap_geo_offset, width=0.05, font.size=2.5, 
               color = NULL,
               hjust = -0.1, colnames_level=colnames(heatmap_3), 
               colnames_angle=0, legend_title="",
               custom_column_labels  = rep("", ncol(heatmap_3))
  ) + 
    scale_fill_manual(values = geo.colours,
                      name = "Geo")+
    theme(legend.position = c(1.1,0.92),
          legend.background = element_rect(fill = "transparent",
                                           colour = NA))
  
  return(p7)
  
  
  
}






# plot_legend<-plot_legend(isotype_map_0.6)
# ggsave("legend.pdf", plot = plot_legend, width = 7.5, height = 7.5, units = "in")

plot_0.6_tree<-plot_tree_heatmaps(tree_file = isotype_map_0.6,
                                  annotation_map_file = annotation_maps,
                                  heatmap_geo_offset= 0.0021)
ggsave("raw_Ct_tree_0.6.pdf", plot = plot_0.6_tree, width = 7.5, height = 7.5, units = "in")

plot_0.7_tree<-plot_tree_heatmaps(tree_file = isotype_map_0.7,
                                  annotation_map_file = annotation_maps,
                                  heatmap_geo_offset= 0.0023)
ggsave("raw_Ct_tree_0.7.pdf", plot = plot_0.7_tree, width = 7.5, height = 7.5, units = "in")

plot_0.8_tree<-plot_tree_heatmaps(tree_file = isotype_map_0.8,
                                  annotation_map_file = annotation_maps)
ggsave("raw_Ct_tree_0.8.pdf", plot = plot_0.8_tree, width = 7.5, height = 7.5, units = "in")

plot_0.9_tree<-plot_tree_heatmaps(tree_file = isotype_map_0.9,
                                  annotation_map_file = annotation_maps)
ggsave("raw_Ct_tree_0.9.pdf", plot = plot_0.9_tree, width = 7.5, height = 7.5, units = "in")

#save image 7.5inches*inches
#rearrange legend positions and titles with Adobe Illustrator











##########  function: plot the equal_angle tree
plot_equal_angle_tree <- function(tree_file, 
                                  annotation_map_file){
  
  annotation_maps<-annotation_map_file
  
  p_ea<-ggtree::ggtree(tree_file, 
                       # aes(col=geo), 
                       layout="equal_angle", size=0.15) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
    geom_tippoint(aes(color=geo), size=1, alpha =0.5, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
    # geom_tiplab(aes(label=label, col=geo), size = 0.5, hjust=-0.5, align=TRUE, linesize=0.1, alpha=0.6)+
    # geom_tiplab(aes(label=NA, col=NA), size=0.8)+
    theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  return(p_ea)
}

annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)
plot_0.9_tree_equal_angle<-plot_equal_angle_tree(isotype_map_0.9_raw, 
                                                 annotation_map_file = annotation_maps)
ggsave("Ct_tree_0.9_unrooted.pdf", plot = plot_0.9_tree_equal_angle, width = 7.5, height = 7.5, units = "in")











##########  function: plot the equal_angle tree for Figure 1 - no legend ###
library(cowplot)
plot_equal_angle_tree_rotated <- function(tree_file, 
                                          annotation_map_file){
  
  annotation_maps<-annotation_map_file
  
  p_ea<-ggtree::ggtree(tree_file, 
                       # aes(col=geo), 
                       layout="equal_angle", size=0.15) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
    geom_tippoint(aes(color=geo), size=1, alpha =0.5, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
    # theme_tree2() +
    # geom_tiplab(aes(label=label, col=geo), size = 0.5, hjust=-0.5, align=TRUE, linesize=0.1, alpha=0.6)+
    # geom_tiplab(aes(label=NA, col=NA), size=0.8)+
    theme(legend.title=element_text(face="bold"), 
          legend.position="none",
          legend.box="horizontal", 
          legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  
  # convert ggplot abject into grob
  g <- ggplotGrob(p_ea)
  rotated_grob <- grid::grobTree(
    g,
    vp = grid::viewport(
      angle = 130, 
      x = unit(0.5, "npc"), 
      y = unit(0.5, "npc"),
      width = unit(1, "npc"),
      height = unit(1, "npc")
    )
  )
  
  # use cowplot to draw
  p_rotated <- cowplot::ggdraw() + draw_grob(rotated_grob)
  
  return(p_rotated)
  
  
  # return(p_ea)
}


plot_0.9_tree_equal_angle_rotated<-plot_equal_angle_tree_rotated(isotype_map_0.9_raw, 
                                                                 annotation_map_file = annotation_maps)

plot_0.9_tree_equal_angle_rotated
saveRDS(plot_0.9_tree_equal_angle_rotated, 
        file = "../../processed_data/assemble_figure_1/Ct_plot_0.9_tree_equal_angle_rotated.rds")

ggsave("Ct_tree_0.9_unrooted_rotated.pdf", plot = plot_0.9_tree_equal_angle_rotated, width = 7.5, height = 7.5, units = "in")
 









