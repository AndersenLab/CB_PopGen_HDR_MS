rm(list = ls())

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(dplyr)
source("../utilities.R")

gtcheck_strain_raw<-read.table("../../processed_data/genetic_similarity_and_admixutre/strain_gtcheck.txt",
                            header = TRUE)
gtcheck_strain<-gtcheck_strain_raw %>%
  dplyr::mutate(concordance = 1-(discordance/sites)) %>%
  dplyr::select(i,j,concordance)

isotype_groups_raw<-read.table("../../data/briggsae_isotypes/isotype_groups.tsv",
                               sep = '\t',header = TRUE)

################### NIC174 ###############
target_isotype<-"NIC174"

strains_target_isotype<-isotype_groups_raw %>% 
  dplyr::filter(isotype %in% target_isotype) %>% 
  dplyr::select(strain) %>% 
  dplyr::pull()

NIC174_strains_output<-c(strains_target_isotype,"SOW22") %>%
  as.data.frame()

gt_matrix <- gtcheck_strain %>%  
  tidyr::pivot_wider(names_from = j, values_from = concordance) %>%
  tibble::column_to_rownames("i") %>%
  as.matrix()

# add new row
new_row <- matrix(NA, nrow = 1, ncol = ncol(gt_matrix))
rownames(new_row) <- colnames(gt_matrix)[ncol(gt_matrix)]
gt_matrix <- rbind(gt_matrix,new_row)

# add new col
new_col <- matrix(NA, ncol = 1, nrow = nrow(gt_matrix))
colnames(new_col) <- rownames(gt_matrix)[1]
gt_matrix <- cbind(new_col,gt_matrix)

# add lower.tri data
gt_matrix[lower.tri(gt_matrix,diag = FALSE)] <- t(gt_matrix)[lower.tri(gt_matrix,diag = FALSE)]

# sort
sorted_rows <- sort(rownames(gt_matrix))
sorted_cols <- sort(colnames(gt_matrix))

gt_matrix_plot_tmp <- gt_matrix[sorted_rows, sorted_cols]
diag(gt_matrix_plot_tmp ) <- 1

gt_matrix_plot<-gt_matrix_plot_tmp %>% 
  as.data.frame() %>% 
  mutate(strain= rownames(gt_matrix_plot_tmp)) %>% 
  filter(strain %in% c(strains_target_isotype,isotype_groups_raw$isotype)) %>% 
  select(all_of(c(strains_target_isotype,isotype_groups_raw$isotype)))%>%
  select(sort(names(.))) %>% 
  as.matrix()

# sort
sorted_rows <- sort(rownames(gt_matrix_plot))
sorted_cols <- sort(colnames(gt_matrix_plot))

########### Annotations #######
### Geo
geo_info_raw<-read.csv("../../processed_data/Geo_info/Cb_indep_strain_info_geo.csv")
geo_info<-geo_info_raw %>%
  dplyr::filter(strain %in% rownames(gt_matrix_plot)) %>% 
  tibble::column_to_rownames(var = "strain") %>% 
  dplyr::select(geo) %>% 
  dplyr::rename(Geo=geo)
geo_info$Geo<-as.factor(geo_info$Geo)

filter_geo.colours <- geo.colours[names(geo.colours) %in% unique(geo_info$Geo)]
col_ordered_geo <- geo_info[sorted_cols, "Geo", drop = FALSE]


#######################################################
##### plot heatmap (add target annotation) ############
#######################################################
col_order <- colnames(gt_matrix_plot)
row_order <- rownames(gt_matrix_plot)
col_target_vec <- factor(ifelse(col_order %in% strains_target_isotype, target_isotype, "Other"),
                         levels = c(target_isotype, "Other"))
row_target_vec <- factor(ifelse(row_order %in% strains_target_isotype, target_isotype, "Other"),
                         levels = c(target_isotype, "Other"))
target_colours <- setNames(c("red", "grey"), c(target_isotype, "Other"))

filter_geo.colours <- geo.colours[names(geo.colours) %in% unique(geo_info$Geo)]
col_ordered_geo <- geo_info[sorted_cols, "Geo", drop = FALSE]

col_ha <- HeatmapAnnotation(
  Geo = col_ordered_geo$Geo,
  Target = col_target_vec,
  col = list(Geo = filter_geo.colours, Target = target_colours),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    Geo = list(title = "Geo",
               title_gp    = gpar(fontsize = 10, fontface = "bold"),
               labels_gp   = gpar(fontsize = 8),
               legend_width  = unit(0.5, "cm"),
               legend_height = unit(2, "cm")
    ),
    Target = list(title = "Target isotype",
                  title_gp  = gpar(fontsize = 10, fontface = "bold"),
                  labels_gp = gpar(fontsize = 8),
                  legend_gp = gpar(fill = c("red","grey")))
  )
)

min_concordance <- min(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)
max_concordance <- max(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)

p_heatmap <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = colorRamp2(
    c(min_concordance,
      0.99,
      max_concordance),
    c("#3182bd", "white", "orange")
  ),
  top_annotation       = col_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows = "average",
  clustering_method_columns = "average", 
  show_row_names       = FALSE,
  show_column_names    = FALSE,
  row_dend_gp          = gpar(lwd = 0.3),
  column_dend_gp       = gpar(lwd = 0.3),
  
  column_title         = paste0(target_isotype," strains and other isotype reference strains"),
  column_title_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_title_rot     = 0
  
)


##### only target isotype
gt_matrix_plot_target<-gt_matrix_plot_tmp %>% 
  as.data.frame() %>% 
  dplyr::mutate(strain= rownames(gt_matrix_plot_tmp)) %>% 
  dplyr::filter(strain %in% strains_target_isotype) %>% 
  dplyr::select(all_of(strains_target_isotype))%>%
  dplyr::select(sort(names(.))) %>% 
  as.matrix()

# sort
sorted_rows_target <- sort(rownames(gt_matrix_plot_target))
sorted_cols_target <- sort(colnames(gt_matrix_plot_target))

min_concordance_target  <- min(gt_matrix_plot_target [gt_matrix_plot_target  != 0], na.rm = TRUE)
max_concordance_target  <- max(gt_matrix_plot_target [gt_matrix_plot_target  != 0], na.rm = TRUE)


#### tmp ### 
tmp_gt_matrix_long <- gt_matrix_plot_target %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>%
  tidyr::pivot_longer(
    cols = -ID,
    names_to = "SNP", 
    values_to = "Value"
  ) %>% 
  dplyr::filter(Value != 1)

############ geo annotation #########
strain_anno_raw<-read.csv("../../processed_data/Geo_info/indep_info_geo_for_each_strain.csv")

## region annotaion
strain_anno<-strain_anno_raw %>% 
  dplyr::filter(strain %in% strains_target_isotype) %>% 
  tibble::column_to_rownames(var = "strain") %>% 
  dplyr::select(strain_geo) 

strain_anno$strain_geo<-as.factor(strain_anno$strain_geo)

# keep colours as ordered geo object
filter_geo.colours_target <- geo.colours[names(geo.colours) %in% unique(strain_anno$strain_geo)]
col_ordered_geo_target <- strain_anno[sorted_cols_target, "strain_geo", drop = FALSE]

## lat, long, geo annotaion
strain_anno <- strain_anno_raw %>%
  dplyr::filter(strain %in% strains_target_isotype) %>%
  tibble::column_to_rownames(var = "strain") %>%
  dplyr::select(strain_geo, lat, long)

# factor region
strain_anno$strain_geo <- as.factor(strain_anno$strain_geo)

col_ordered_geo_target  <- strain_anno[sorted_cols_target, "strain_geo", drop = FALSE]
col_ordered_lat_target  <- as.numeric(strain_anno[sorted_cols_target, "lat"])
col_ordered_long_target <- as.numeric(strain_anno[sorted_cols_target, "long"])

# colours 
filter_geo.colours_target <- geo.colours[names(geo.colours) %in% unique(col_ordered_geo_target$strain_geo)]

lat_col_fun <- colorRamp2(
  c(38, mean(c(38, 52),na.rm = TRUE), 52),
  c("whitesmoke", "#6baed6", "#08306b")
)

long_col_fun <- colorRamp2(
  c(-2, mean(c(-2,25), na.rm = TRUE), 25),
  c("whitesmoke", "#fa9752", "#67000d")
)


# define column annotation: Geo + Lat + Long
col_ha_target <- HeatmapAnnotation(
  Geo  = col_ordered_geo_target$strain_geo,
  Lat  = col_ordered_lat_target,
  Long = col_ordered_long_target,
  col = list(
    Geo  = filter_geo.colours_target,
    Lat  = lat_col_fun,
    Long = long_col_fun
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(
    Geo = list(
      title = "Geo",
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_width  = unit(0.6, "cm"),
      legend_height = unit(2, "cm")
    ),
    Lat = list(
      title = "Latitude",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
      # ,at = pretty(lat_rng, n = 3)
    ),
    Long = list(
      title = "Longitude",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
      # ,
      # at = pretty(long_rng, n = 3)
    )
  )
)


#### threshold
cc_iso <- gtcheck_strain %>% 
  dplyr::filter((i %in% strains_target_isotype) & (j %in% strains_target_isotype))

ggplot()+ geom_histogram(data=cc_iso,aes(x=concordance))


library(grid)
p_heatmap <- Heatmap(
  matrix               = gt_matrix_plot_target,
  name                 = "Concordance",
  col                  = colorRamp2(
    c(min_concordance_target,
      0.9999,
      max_concordance_target),
    c( "#3182bd", "white", "orange")
  ),
  top_annotation       = col_ha_target,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3),
  
  
  column_names_rot     = 90,
  column_names_gp      = gpar(fontsize = 1.5),
  column_names_side    = "top",
  column_names_max_height = unit(0.15, "cm"),
  
  
  column_title         = paste0(target_isotype," strains"),
  column_title_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_title_rot     = 0

)


##### only target isotype ##
##### If add SOW22 ###############

strains_target_isotype<-c("SOW22",strains_target_isotype)

gt_matrix_plot_target<-gt_matrix_plot_tmp %>% 
  as.data.frame() %>% 
  mutate(strain= rownames(gt_matrix_plot_tmp)) %>% 
  filter(strain %in% strains_target_isotype) %>% 
  select(all_of(strains_target_isotype))%>%
  select(sort(names(.))) %>% 
  as.matrix()

# sort
sorted_rows_target <- sort(rownames(gt_matrix_plot_target))
sorted_cols_target <- sort(colnames(gt_matrix_plot_target))

min_concordance_target  <- min(gt_matrix_plot_target [gt_matrix_plot_target  != 0], na.rm = TRUE)
max_concordance_target  <- max(gt_matrix_plot_target [gt_matrix_plot_target  != 0], na.rm = TRUE)


############ geo annotation #########
strain_anno_raw<-read.csv("../../processed_data/Geo_info/indep_info_geo_for_each_strain.csv")

## region annotaion
strain_anno<-strain_anno_raw %>% 
  dplyr::filter(strain %in% strains_target_isotype) %>% 
  tibble::column_to_rownames(var = "strain") %>% 
  dplyr::select(strain_geo) 

strain_anno$strain_geo<-as.factor(strain_anno$strain_geo)

# keep colours as ordered geo object
filter_geo.colours_target <- geo.colours[names(geo.colours) %in% unique(strain_anno$strain_geo)]
col_ordered_geo_target <- strain_anno[sorted_cols_target, "strain_geo", drop = FALSE]

## lat, long, geo annotaion
strain_anno <- strain_anno_raw %>%
  filter(strain %in% strains_target_isotype) %>%
  tibble::column_to_rownames(var = "strain") %>%
  dplyr::select(strain_geo, lat, long)

# factor region
strain_anno$strain_geo <- as.factor(strain_anno$strain_geo)
col_ordered_geo_target  <- strain_anno[sorted_cols_target, "strain_geo", drop = FALSE]
col_ordered_lat_target  <- as.numeric(strain_anno[sorted_cols_target, "lat"])
col_ordered_long_target <- as.numeric(strain_anno[sorted_cols_target, "long"])

# colours 
filter_geo.colours_target <- geo.colours[names(geo.colours) %in% unique(col_ordered_geo_target$strain_geo)]

lat_col_fun <- colorRamp2(
  c(38, mean(c(38, 52),na.rm = TRUE), 52),
  c("whitesmoke", "#6baed6", "#08306b")
)

long_col_fun <- colorRamp2(
  c(-2, mean(c(-2,25), na.rm = TRUE), 25),
  c("whitesmoke", "#fa9752", "#67000d")
)

# define column annotation: Geo + Lat + Long
col_ha_target <- HeatmapAnnotation(
  Geo  = col_ordered_geo_target$strain_geo,
  Lat  = col_ordered_lat_target,
  Long = col_ordered_long_target,
  col = list(
    Geo  = filter_geo.colours_target,
    Lat  = lat_col_fun,
    Long = long_col_fun
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(
    Geo = list(
      title = "Geo",
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_width  = unit(0.6, "cm"),
      legend_height = unit(2, "cm")
    ),
    Lat = list(
      title = "Latitude",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
      # ,at = pretty(lat_rng, n = 3)
    ),
    Long = list(
      title = "Longitude",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
    )
  )
)

#### threshold
cc_iso <- gtcheck_strain %>% 
  dplyr::filter((i %in% strains_target_isotype) & (j %in% strains_target_isotype))

# ggplot()+ geom_histogram(data=cc_iso,aes(x=concordance))

p_heatmap <- Heatmap(
  matrix               = gt_matrix_plot_target,
  name                 = "Concordance",
  col                  = colorRamp2(
    c(min_concordance_target,
      0.9999,
      max_concordance_target),
    c( "#3182bd", "white", "orange")
  ),
  top_annotation       = col_ha_target,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3),
  column_names_rot     = 90,
  column_names_gp      = gpar(fontsize = 1.5),
  column_names_side    = "top",
  column_names_max_height = unit(0.15, "cm"),
  
  
  column_title         = paste0(target_isotype," strains with the closet"),
  column_title_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_title_rot     = 0
  
)

p_heatmap

dev.off()
pdf(paste0("../../figures/EDF2_Cb_",target_isotype,"_concordance_heatmap_with_the_closest.pdf"), width = 7, height = 6)


draw(p_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()


