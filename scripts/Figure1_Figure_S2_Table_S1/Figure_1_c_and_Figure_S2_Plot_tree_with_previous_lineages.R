
rm(list=ls())

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ggnewscale)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tibble)
library(readxl)
library(ape)
library(phytools)



library(RColorBrewer)


source("../utilities.R")



display.brewer.all()
display.brewer.pal(11,"RdYlBu")
brewer.pal(11,"RdYlBu")
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")
display.brewer.pal(8,"Set2")
brewer.pal(8,"Set2")
display.brewer.pal(10,"Set3")
brewer.pal(10,"Set3")
display.brewer.pal(8,"Dark2")
brewer.pal(8,"Dark2")




# 1. load tree
isotype_map_0.9<-treeio::read.tree("../../processed_data/test_GTR_LD09/phy_file_LD_0.9.phy.contree")

###### 2. use mid point to Mark the root of the tree ###
isotype_map_0.9<-phytools::midpoint_root(isotype_map_0.9)



### load geo
annotation_maps_raw<- read.csv("../../processed_data/geo_info/Cb_indep_isotype_info_geo.csv", header = TRUE)



####################################################
#######             load raw data           ############
raw_data<-read.csv("../../data/20250626_c_briggsae_strain_data.csv")
indep_isotype_info<- raw_data

sample_list<- read.table("../../processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt")
nrow(sample_list)
#719

WI_isotype_info<-raw_data %>% 
  filter(strain %in% sample_list$V1) %>% 
  rename(lat=latitude) %>% 
  rename(long=longitude) %>% 
  mutate(isotype=strain) %>% 
  select(isotype,lat,long)



############################################################
######## calculate average temperature of each sample ######

### load Hawaii in situ temperature data
library(readxl)
raw_excel<-readxl::read_excel("../../data/C. briggsae WI strain info_20250515.xlsx")

sample_list<- read.table("../../processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt")
nrow(sample_list)
#719

WI_isotype_info_temperature<-raw_excel %>% 
  filter(strain %in% sample_list$V1) %>% 
  rename(lat=latitude) %>% 
  rename(long=longitude) %>% 
  mutate(isotype=strain) %>% 
  select(isotype,lat,long,ambient_temp) %>% 
  mutate(ambient_temp = na_if(ambient_temp, "NA"))  ### There is a lot of manually typed "NA" characters in the raw data sheet

WI_isotype_info_temperature<-WI_isotype_info_temperature

### load NOAA temperature results
NOAA_3_month_raw<-read.delim("../../processed_data/NOAA/Cb_NOAA_from_Mike/Cb_3month_NOAA/phenotypes.tsv",
                             sep = '\t',
                             header = TRUE)

NOAA_3_month<-NOAA_3_month_raw %>%
  select(isotype,mean_mean.TEMP) %>% 
  rename(NOAA_temp=mean_mean.TEMP)


average_temperature<-WI_isotype_info_temperature %>% 
  left_join(NOAA_3_month, by= c("isotype") ) %>%
  mutate(ambient_temp = ifelse(is.na(ambient_temp),
                               NOAA_temp,ambient_temp)) %>%
  rename(average_temp=ambient_temp) %>%
  select(-NOAA_temp) %>%
  mutate(average_temp = ifelse(average_temp=="na", NA ,average_temp)) %>% 
  mutate(average_temp = as.numeric(average_temp)) %>%
  as.data.frame()

write.table(average_temperature,
            "isotypes_average_temperature.tsv",
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)


######### End - calculate average temperature of each sample ##########
#########################################################################










#### Thomas paper Lineages annotations
thomas_paper_tree_annotation_raw<-read.table("../../tables/Thomas_paper_tree_annotation/thomas_paper_tree_annotation.txt",
                 header = TRUE)
thomas_paper_tree_annotation<-thomas_paper_tree_annotation_raw %>% 
  rename(Lineages=Phylogeographic_clade) %>% 
  rename(Lineages_color = color) %>% 
  select(isotype, Lineages, Lineages_color) %>% 
  na.omit() %>% 
  unique() %>% 
  filter(isotype != "NIC174")
  # Remove NIC174, Because Both JU1376 and JU516 were assigned in the same isotype NIC174. 
  # But they were defined to two Lineages, temperate and tropical, as defined by Thomas et al. 2015.


annotation_Thomas_color<- thomas_paper_tree_annotation %>%
  distinct(Lineages, Lineages_color) %>%
  tibble::deframe()







### Marie-anne and Richard paper Lineages annotations
gene_segments_trees_raw<-read_excel("../../data/From_Marie_anne_2013_and_Richard_2011_paper/all_strains_manually_generated.xlsx")

gene_segments_trees<-gene_segments_trees_raw %>% 
  rename(Lineage = Clade) %>% 
  mutate(Lineage = ifelse(Lineage == "Tropical I", "Tropical", Lineage)) %>% 
  mutate(Lineage = ifelse(Lineage == "Temperate II", "Temperate", Lineage)) %>% 
  mutate(Strain = ifelse(Strain == "VX0033", "VX33", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "VX0034", "VX34", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "GXW0021", "GXW21", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "GXW0022", "GXW22", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "GXW0023", "GXW23", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "GXW0024", "GXW24", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "GXW0025", "GXW25", Strain)) %>% 
  mutate(Strain = ifelse(Strain == "GXW0026", "GXW26", Strain)) %>% 
  rename(strain=Strain)

gene_segments_trees_annotation_tmp<-indep_isotype_info %>% 
  left_join(gene_segments_trees, by = "strain") %>% 
  select(strain, isotype, Lineage, From_which_paper) %>% 
  na.omit()


setdiff(gene_segments_trees$strain,
        gene_segments_trees_annotation_tmp$strain)
# [1] "NIC85"  "NIC90"  "NIC102" "QX1799" "QX1805" "JU376" 
# [7] "VX33"   "ED3079" "ED3084" "ED3099" "ED3100"



## Annotation
gene_segments_trees_annotation<-gene_segments_trees_annotation_tmp %>% 
  select(-strain,-From_which_paper) %>% 
  unique()


strain_in_same_isotype_but_diff_Lineage_MF_paper<-gene_segments_trees_annotation %>% 
  group_by(isotype) %>%
  summarise(count = n()) %>% 
  filter(count >1)

strain_in_same_isotype_but_diff_Lineage_MF_paper
# # A tibble: 2 × 2
# isotype count
# <chr>   <int>
#   1 JU1344      3
# 2 NIC174      2


## now remove them 
gene_segments_trees_annotation<-gene_segments_trees_annotation %>% 
  filter(!(isotype %in% c("JU1344","NIC174")))




### update annotation file, add Thomas paper and Marie-anne paper Lineages
annotation_maps<-annotation_maps_raw %>% 
  left_join(thomas_paper_tree_annotation, by = "isotype") %>% 
  select(-Lineages_color) %>% 
  left_join(gene_segments_trees_annotation, by = "isotype")




############ plot the tree with three heatmaps #########

# function: plot the tree
plot_tree_three_heatmap <- function(tree_file, 
                                  annotation_map_file,
                                  annotation_temperature_file,
                                  annotation_Thomas_file,
                                  annotation_MAF_file,
                                  tree_layout="circular",
                                  add_tiplab=TRUE,
                                  heatmap_lat_offset = 0.0018,
                                  heatmap_temp_offset = 0.0007,
                                  # heatmap_geo_offset = 0.00055,
                                  heatmap_geo_offset = 0.0029,
                                  heatmap_thomas_offset = 0.0040,
                                  heatmap_MAF_offset = 0.0051){
  
  data_tree_file<-ggtree::fortify(tree_file)
  annotation_maps<- annotation_map_file
  annotation_temperature <- annotation_temperature_file
  annotation_Thomas<-annotation_Thomas_file
  annotation_MAF<-annotation_MAF_file
  
  heatmap<- annotation_maps%>%
    select(lat)
  row.names(heatmap)<-annotation_maps$isotype
  heatmap$lat<-abs(heatmap$lat)
  
  heatmap_2<- annotation_temperature%>%
    select(average_temp)
  row.names(heatmap_2)<-annotation_temperature$isotype
  
  heatmap_3<- annotation_maps%>%
    select(geo)
  row.names(heatmap_3)<-annotation_maps$isotype
  
  
  
  heatmap_4 <- annotation_Thomas %>%
    distinct(isotype, Lineages) %>%  
    column_to_rownames(var = "isotype") %>% 
    select(Lineages)

  annotation_Thomas_color<- annotation_Thomas %>%
    distinct(Lineages, Lineages_color) %>%
    tibble::deframe()
  
  
  heatmap_5 <- annotation_MAF %>%
    column_to_rownames(var = "isotype") %>% 
    select(Lineage)
  

  
  # plot tree
  p0<-ggtree::ggtree(data_tree_file, aes(col=geo), layout=tree_layout, size=0.15) %<+% annotation_maps +
    geom_tippoint(aes(color=geo), size=0.1, alpha =0.6, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
 
    theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  if (add_tiplab==TRUE) {
    p0<- p0+geom_tiplab(aes(label=label),
                color = "grey50",
                linetype = "dashed",
                size = 0.5, hjust = -0.5,
                align = TRUE, linesize = 0.05,
                alpha = 0.6)
  }
  
  
  p1<- ggtree::open_tree(p0, 90) %>% ggtree::rotate_tree(90) 
  
  
  p2<-p1+new_scale_fill()+new_scale_color()
  p3<-gheatmap(p2, heatmap, offset = heatmap_lat_offset, width=0.05, font.size=2.5, 
               color = NULL, hjust = -0.1, colnames_level=colnames(heatmap), 
               colnames_angle=0, 
               legend_title=" Absolute latitude",
               custom_column_labels  = rep(" Abs lat", ncol(heatmap))
  ) + 
    scale_fill_gradient2(
      low = "#D73027", mid = "#FFFFBF",high = "#4575B4",
      midpoint = mean(c(max(heatmap$lat, na.rm = TRUE), min(heatmap$lat, na.rm = TRUE))),
      na.value = "white",  
      name = "Absolute latitude"
    )

  p4<-p3+new_scale_fill()+new_scale_color()
 
    p6<-p4+new_scale_fill()+new_scale_color()
    p7<-gheatmap(p6, heatmap_3, offset = heatmap_geo_offset, width=0.05, font.size=2.5, 
                 color = NULL,
                 hjust = -0.1, colnames_level=colnames(heatmap_3), 
                 colnames_angle=0, legend_title="",
                 custom_column_labels  = rep(" Geo", ncol(heatmap_3))
    ) + 
      scale_fill_manual(values = geo.colours,
                         name = "Geo")+
      theme(legend.position = c(0.99,0.72),
            legend.background = element_rect(fill = "transparent",
                                             colour = NA))
    
    
    p8<-p7+new_scale_fill()+new_scale_color()
    p9<-gheatmap(p8, heatmap_4, offset = heatmap_thomas_offset, width=0.05, font.size=2.5, 
                 color = NULL,
                 hjust = -0.1, colnames_level=colnames(heatmap_4), 
                 colnames_angle=0, legend_title="",
                 custom_column_labels  = rep("2015 Lineage", ncol(heatmap_4))
    ) + 
      scale_fill_manual(values = annotation_Thomas_color,
                        name = "Lineage",
                        na.value = "white")+
      theme(legend.position = c(0.99,0.72),
            legend.background = element_rect(fill = "transparent",
                                             colour = NA))
    
    
    p10<-p9+new_scale_fill()+new_scale_color()
    p11<-gheatmap(p10, heatmap_5, offset = heatmap_MAF_offset, width=0.05, font.size=2.5, 
                 color = NULL,
                 hjust = -0.1, colnames_level=colnames(heatmap_5), 
                 colnames_angle=0, legend_title="",
                 custom_column_labels  = rep("2013 Lineage", ncol(heatmap_5))
    ) + 
      scale_fill_manual(values = annotation_Thomas_color,
                        name = NULL,
                        na.value = "white")+
      theme(legend.position = c(0.99,0.72),
            legend.background = element_rect(fill = "transparent",
                                             colour = NA))
    
    p12<-p11+new_scale_fill()+new_scale_color()
    p13<-gheatmap(p12, heatmap_2, offset = heatmap_temp_offset, width=0.05, font.size=2.5, 
                 color = NULL, hjust = -0.1, colnames_level=colnames(heatmap_2), 
                 colnames_angle=0, legend_title=" Mean temperature °C",
                 custom_column_labels  = rep("Mean temp", ncol(heatmap_2))
    ) + 
      scale_fill_gradient2(
        low = "#4575B4", mid = "#FFFFBF",high = "#D73027",
        midpoint = mean(c(max(heatmap_2$average_temp, na.rm = TRUE), min(heatmap_2$average_temp, na.rm = TRUE))),
        na.value = "white",  # don't show NA in the heatmap
        name = "Mean temperature °C"
      )+
      theme(legend.position = c(1.2,0.72),
            legend.background = element_rect(fill = "transparent",
                                             colour = NA))
    

    
}





# plot_0.6_tree_2heatmap<-plot_tree_three_heatmap(isotype_map_0.6, annotation_map_file = annotation_maps,annotation_temperature_file = average_temperature)
# ggsave("raw_Cb_tree_0.6_2heatmap.pdf", plot = plot_0.6_tree_2heatmap, width = 7.5, height = 7.5, units = "in")

plot_0.9_tree_2heatmap_main<-plot_tree_three_heatmap(isotype_map_0.9, 
                                                annotation_map_file = annotation_maps,
                                                annotation_temperature_file = average_temperature,
                                                annotation_Thomas_file = thomas_paper_tree_annotation,
                                                annotation_MAF_file = gene_segments_trees_annotation,
                                                add_tiplab = FALSE)



ggsave("raw_Cb_tree_0.9_3heatmap_main.pdf", plot = plot_0.9_tree_2heatmap_main, width = 7, height = 7, units = "in")






############ Adding label - plot the tree with heatmaps #########

plot_0.9_tree_2heatmap_supp<-plot_tree_three_heatmap(isotype_map_0.9, 
                                               annotation_map_file = annotation_maps,
                                               annotation_temperature_file = average_temperature,
                                               annotation_Thomas_file = thomas_paper_tree_annotation,
                                               annotation_MAF_file = gene_segments_trees_annotation,
                                               # add_tiplab = FALSE
)
plot_0.9_tree_2heatmap_supp
ggsave("raw_Cb_tree_0.9_3heatmap_supp.pdf", 
       plot = plot_0.9_tree_2heatmap_supp, 
       width = 7, height = 7, units = "in")





















##########  function: plot the equal_angle tree for Figure 1 - no legend ###
plot_equal_angle_tree_rotated <- function(tree_file, 
                                          annotation_map_file){
  
  annotation_maps<-annotation_map_file
  
  p_ea<-ggtree::ggtree(tree_file, 
                       # aes(col=geo), 
                       layout="equal_angle", size=0.05) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
    geom_tippoint(aes(color=geo), size=1, alpha =0.5, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
 
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
      angle = -55, 
      x = unit(0.5, "npc"), 
      y = unit(0.5, "npc"),
      width = unit(1, "npc"),
      height = unit(1, "npc")
    )
  )
  
  # use cowplot to draw
  p_rotated <- ggdraw() + draw_grob(rotated_grob)
  
  return(p_rotated)
  
  
  # return(p_ea)
}


plot_0.9_tree_equal_angle_rotated<-plot_equal_angle_tree_rotated(isotype_map_0.9, 
                                                                 annotation_map_file = annotation_maps)

plot_0.9_tree_equal_angle_rotated
saveRDS(plot_0.9_tree_equal_angle_rotated, 
        file = "../../processed_data/assemble_figure_1/Cb_plot_0.9_tree_equal_angle_rotated.rds")

ggsave("Cb_tree_0.9_unrooted_rotated.pdf", 
       plot = plot_0.9_tree_equal_angle_rotated, 
       width = 7, height = 7, units = "in")

