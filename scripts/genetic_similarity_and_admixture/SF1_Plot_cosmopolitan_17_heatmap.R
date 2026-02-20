rm(list = ls())

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(dplyr)
library(gridExtra)
library(grid)
source("../utilities.R")

gtcheck_isotype_raw<-read.table("../../processed_data/heatmap_hard_filtered/gtcheck.tsv",
                               header = TRUE)
gtcheck_isotype<-gtcheck_isotype_raw %>%
  dplyr::mutate(concordance = 1-(discordance/sites)) %>%
  dplyr::select(i,j,concordance)

################### 17 isotypes ###############
target_18_isotype<-read.table("../../processed_data/Geo_info/Cb_cosmopolitan_isotype.txt") %>% 
  dplyr::select(V1) %>% 
  dplyr::pull()

gt_matrix <-gtcheck_isotype %>%
  tidyr::pivot_wider(names_from = j, values_from = concordance) %>%
  tibble::column_to_rownames("i") %>%
  as.matrix()

# add new row
new_row <- matrix(NA, nrow = 1, ncol = ncol(gt_matrix))
rownames(new_row) <- colnames(gt_matrix)[1]
gt_matrix <- rbind(new_row,gt_matrix)

# add new col
new_col <- matrix(NA, ncol = 1, nrow = nrow(gt_matrix))
colnames(new_col) <- rownames(gt_matrix)[nrow(gt_matrix)]
gt_matrix <- cbind(gt_matrix,new_col)

# add upper.tri data
gt_matrix[upper.tri(gt_matrix,diag = FALSE)] <- t(gt_matrix)[upper.tri(gt_matrix,diag = FALSE)]

# sort
sorted_rows <- sort(rownames(gt_matrix))
sorted_cols <- sort(colnames(gt_matrix))

gt_matrix_plot_tmp <- gt_matrix[sorted_rows, sorted_cols]
diag(gt_matrix_plot_tmp ) <- 1

gt_matrix_plot<-gt_matrix_plot_tmp %>% 
  as.data.frame() %>% 
  dplyr::mutate(isotype = rownames(gt_matrix_plot_tmp)) %>% 
  dplyr::filter(isotype %in% target_18_isotype) %>% 
  dplyr::select(all_of(target_18_isotype))%>%
  dplyr::select(sort(names(.))) %>% 
  as.matrix()

# sort
sorted_rows <- sort(rownames(gt_matrix_plot))
sorted_cols <- sort(colnames(gt_matrix_plot))

min_concordance  <- min(gt_matrix_plot [gt_matrix_plot  != 0], na.rm = TRUE)
max_concordance  <- max(gt_matrix_plot [gt_matrix_plot  != 0], na.rm = TRUE)

##############################
##### no-annotation ######
##############################

phylo_vals <- as.vector(gt_matrix_plot)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]  # remove NAs

breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)

phylo_col_fun <- circlize::colorRamp2(
  breaks,
  c("#4575B4", "#87C6C2", "#FFFFE0", "#F4D166", "#D73027")
)

p_heatmap_1 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Genetic\nsimilarity",
  col                  = phylo_col_fun,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows    = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows      = "average",
  clustering_method_columns   = "average",
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 6),
  column_names_side    = "top",
  row_dend_gp          = gpar(lwd = 0.3),
  column_dend_gp       = gpar(lwd = 0.3)
)

##############################
##### lineage-annotation ######
##############################
lin_colors <- c(TI="#ff77ab", TC="#ff0000", Tropical="#ff0000",TT="#ff0000",NWD="#00ff00",#WD="#d0fe02",
                Montreal="#ff037e", Hubei="#16537e", ID="#6a329f", TD1="#f3c588", TD2="#ffdd02",
                TD3="#00f4c2", KD = "#8b3700", Temperate="#0000ff", TA="#9fc5e8", AD="#ff8200", TH="#aeb400",FM="#000000",
                Quebec="#ff037e")

lineage_raw<-readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv")

lineage_18<-lineage_raw %>% 
  dplyr::filter(isotype %in% target_18_isotype) %>% 
  dplyr::select(isotype, Lineage)

present_lineages <- unique(lineage_18$Lineage)
lin_colors_filtered <- lin_colors[names(lin_colors) %in% present_lineages]

lineage_df <- lineage_18 %>%
  dplyr::filter(isotype %in% colnames(gt_matrix_plot)) %>%
  dplyr::arrange(match(isotype, colnames(gt_matrix_plot)))

col_ha <- HeatmapAnnotation(
  group = lineage_df$Lineage,
  col = list(group = lin_colors_filtered), 
  annotation_name_side = "left", 
  annotation_legend_param = list(
    group = list(
      title = "Relatedness\ngroup", 
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 6)
    )
  ),
  show_annotation_name = TRUE,
  show_legend = TRUE
)

phylo_vals <- as.vector(gt_matrix_plot)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]

breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)

phylo_col_fun <- circlize::colorRamp2(
  breaks,
  c("#4575B4", "#87C6C2", "#FFFFE0", "#F4D166", "#D73027")
)

p_heatmap_2 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Genetic\nsimilarity",
  col                  = phylo_col_fun,
  top_annotation       = col_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows    = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows      = "average",
  clustering_method_columns   = "average",
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 6) 
  ),
  column_names_gp      = gpar(fontsize = 6),
  column_names_side    = "top",
  row_dend_gp          = gpar(lwd = 0.3),
  column_dend_gp       = gpar(lwd = 0.3)
)



######## assemble Figure S2 #########
load("../../processed_data/genetic_similarity_and_admixutre/Cosmo_18_maps.RData")

g1 <- ggplotGrob(combined_plot)
g2 <- grid.grabExpr(draw(p_heatmap_1))

dev.off()
pdf("../../figures/raw_SF1_cosmopolitan_isotypes.pdf", width = 7, height = 8)

grid.arrange(g1, g2, ncol = 1, heights = c(1, 1))

grid.text("a", 
          x = unit(0, "npc") + unit(0.5, "mm"), 
          y = unit(1, "npc") - unit(2, "mm"), 
          just = c("left","top"),
          gp = gpar(fontface = "bold",fontsize = 14))
grid.text("b", 
          x = unit(0, "npc") + unit(0.5, "mm"), 
          y = unit(0.45, "npc") - unit(2, "mm"), 
          just = c("left","top"),
          gp = gpar(fontface = "bold",fontsize = 14))

dev.off()

