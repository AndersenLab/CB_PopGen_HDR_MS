rm(list=ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(dendextend)
library(ggpubr)
library(cowplot)
library(patchwork)
library(ggridges)
library(readr)

source("../utilities.R")

geo_info_raw<-read.csv("../../processed_data/geo_info/Cb_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw 

lineage_raw<-read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv") %>% 
  dplyr::filter(!(isotype %in% c("MY681", "ECA1146", "JU356", "ECA1503")))

plot_PCA<-function(PCA_input, tracy_for_plot_input, x_axis, y_axis){

  label_x<-tracy_for_plot_input %>%
    dplyr::filter(N == as.numeric(sub("PC", "", x_axis))) %>%
    dplyr::select(VarExp) %>%
    pull(VarExp) %>%
    `*`(100) %>%
    round(digits = 2) %>%
    as.character()

  label_y<-tracy_for_plot_input %>%
    dplyr::filter(N == as.numeric(sub("PC", "", y_axis))) %>%
    dplyr::select(VarExp) %>%
    pull(VarExp) %>%
    `*`(100) %>%
    round(digits = 2) %>%
    as.character()

  p1_1<-ggplot2::ggplot(PCA_input)+
    geom_point(shape=16, alpha=0.8, size=1.2, aes(x=.data[[x_axis]],y=.data[[y_axis]], color=geo))+
    scale_color_manual(values = geo.colours, name = "geo") +
    theme_bw() +
    labs(x=paste(x_axis," (",label_x,"%)",sep = ""),
         y=paste(y_axis," (",label_y,"%)",sep = ""))+
    theme(axis.title = element_text(size=6, face = "bold",color = "black"),
          axis.text = element_text(size=5, color = "black"),
          panel.grid = element_blank())+
    theme(legend.position = "none")
 
  p1_1

  return(p1_1)
}

TracyWidom_file<-"../../processed_data/PCA_more_iterations/PCA_5_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_file<-"../../processed_data/PCA_more_iterations/PCA_5_iterations/eigenstrat_outliers_removed.evac"

plot_remove_out_liner<-function(TracyWidom_file,
                                eigenstrat_file,
                                n_iterations){

# calculate % var explained.
tracy_for_plot_rm_out_liner <- data.table::fread(TracyWidom_file) %>%
  dplyr::mutate(sum = sum(eigenvalue),
                VarExp = eigenvalue/sum,
                sigEV = ifelse(`p-value` < 0.05, T, F)) %>%
  dplyr::group_by(sigEV) %>%
  dplyr::mutate(sigVarExp = sum(VarExp)) %>%
  dplyr::ungroup()

# make priciple components df with annotated collections 0.9
cdf9no_rm_out_liner <- data.table::fread(eigenstrat_file, skip = 1)

cdf9_2no_rm_out_liner <- cdf9no_rm_out_liner %>%
  dplyr::select(isotype = V1, PC1=V2, PC2=V3, PC3=V4, PC4=V5, PC5=V6, PC6 = V7) #  six significant eigenvectors wihtout outlier removal

cdf9_3no_rm_out_liner <- cdf9_2no_rm_out_liner %>% dplyr::select(-isotype)

# add back location data nd use these dfs to plot PC1 by PC2
cdf9_7no_rm_out_liner <- cdf9_2no_rm_out_liner %>%
  dplyr::left_join(geo_info,by=c("isotype")) 

pca_TAC_ld0.9_rm_out_liner <- cdf9_7no_rm_out_liner
n_of_isotype<-nrow(pca_TAC_ld0.9_rm_out_liner)

p_PC1_PC2_rm_out_liner<- plot_PCA(pca_TAC_ld0.9_rm_out_liner,
                                  tracy_for_plot_input = tracy_for_plot_rm_out_liner,
                                  x_axis="PC1",y_axis="PC2")

p_PC3_PC4_rm_out_liner<- plot_PCA(PCA_input=pca_TAC_ld0.9_rm_out_liner,
                                  tracy_for_plot_input = tracy_for_plot_rm_out_liner,
                                  x_axis="PC3",y_axis="PC4")

# generate legend
legend_plot_rm_out_liner <- ggplot(pca_TAC_ld0.9_rm_out_liner)+
  geom_point(aes(x = PC1, y = PC2, color = geo), alpha=0.8, size=1.2, shape = 16) +
  scale_color_manual(values = geo.colours, name = "Geo") +
  theme_bw() +
  guides(color = guide_legend(
    nrow = 2,
    byrow = TRUE,
    keywidth = unit(0.8, "lines"),
    keyheight = unit(0.8, "lines")
  )) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.spacing = unit(0.05, "cm"),  
    legend.spacing.x = unit(0.08, "cm"), 
    legend.spacing.y = unit(0.01, "cm"),   
    legend.key.size = unit(0.8, "lines"), 
    legend.text = element_text(size = 7),  
    legend.title = element_text(size = 8),  
    legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
  )

legend_density_rm_out_liner <-cowplot::get_plot_component(legend_plot_rm_out_liner, "guide-box-bottom")

Cb_rm_out_liner_PCA_plot_PC1_4_with_density <- (
  p_PC1_PC2_rm_out_liner | p_PC3_PC4_rm_out_liner  
) 

Cb_rm_out_liner_PCA_plot_PC1_4_with_density<-Cb_rm_out_liner_PCA_plot_PC1_4_with_density +
  plot_annotation(
    title = paste0("PCA after ",
                   n_iterations,
                   " iterations of outlier removal  -  n_isotypes = ",
                   n_of_isotype)
  ) &
  theme(
    plot.title = element_text(size = 7, face = "bold", 
                              hjust = 0.5, vjust = -4),
  )

return(list(p_iteration=Cb_rm_out_liner_PCA_plot_PC1_4_with_density,
            common_legend=legend_density_rm_out_liner,
            n_of_isotype=n_of_isotype,
            pca_TAC_ld0.9_rm_out_liner=pca_TAC_ld0.9_rm_out_liner
       )
)
}

#### 5 iterations
TracyWidom_5<-"../../processed_data/PCA_more_iterations/PCA_5_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_5<-"../../processed_data/PCA_more_iterations/PCA_5_iterations/eigenstrat_outliers_removed.evac"

iteration_5_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_5,
                                        eigenstrat_file = eigenstrat_5,
                                        n_iterations = "5")

iter_5_df<-iteration_5_results$pca_TAC_ld0.9_rm_out_liner %>% 
  dplyr::left_join(lineage_raw %>% select(isotype,Lineage), by = c("isotype"))

PC12_up<-iter_5_df %>% 
  dplyr::filter(PC2 >0.15)

PC12_down_big<-iter_5_df %>% 
  dplyr::filter(PC2 <0.05) %>% 
  dplyr::filter(PC1 <0.2)%>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())
  
PC34_3over_0<-iter_5_df %>% 
  dplyr::filter(PC3 >0) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC34_left<-iter_5_df %>% 
  dplyr::filter(PC3 < -0.05) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#removed isotypes in the first 5 iterations
diff_iso<-setdiff(lineage_raw$isotype,iter_5_df$isotype)

removed_isotypes_df<-lineage_raw %>% 
  dplyr::filter(isotype %in% diff_iso)

removed_isotypes_df_summary<-removed_isotypes_df %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#### 10 iterations
TracyWidom_10<-"../../processed_data/PCA_more_iterations/PCA_10_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_10<-"../../processed_data/PCA_more_iterations/PCA_10_iterations/eigenstrat_outliers_removed.evac"

iteration_10_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_10,
                                        eigenstrat_file = eigenstrat_10,
                                        n_iterations = "10")

iter_10_df<-iteration_10_results$pca_TAC_ld0.9_rm_out_liner %>% 
  dplyr::left_join(lineage_raw %>% select(isotype,Lineage), by = c("isotype"))

PC12_up<-iter_10_df %>% 
  dplyr::filter(PC2 >0.04) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_right<-iter_10_df %>% 
  dplyr::filter(PC1 >0.15) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_mid<-iter_10_df %>% 
  dplyr::filter(PC1 <0.15 & PC1 >0.05) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#### 15 iterations
TracyWidom_15<-"../../processed_data/PCA_more_iterations/PCA_15_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_15<-"../../processed_data/PCA_more_iterations/PCA_15_iterations/eigenstrat_outliers_removed.evac"

iteration_15_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_15,
                                         eigenstrat_file = eigenstrat_15,
                                         n_iterations = "15")

iter_15_df<-iteration_15_results$pca_TAC_ld0.9_rm_out_liner %>% 
  dplyr::left_join(lineage_raw %>% select(isotype,Lineage), by = c("isotype"))

PC12_up<-iter_15_df %>% 
  dplyr::filter(PC2>0.04) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_left<-iter_15_df %>% 
  dplyr::filter(PC2<0.04) %>% 
  dplyr::filter(PC1<0.05) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())


PC12_down_right<-iter_15_df %>% 
  dplyr::filter(PC1>0.2) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#### 20 iterations
TracyWidom_20<-"../../processed_data/PCA_more_iterations/PCA_20_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_20<-"../../processed_data/PCA_more_iterations/PCA_20_iterations/eigenstrat_outliers_removed.evac"

iteration_20_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_20,
                                         eigenstrat_file = eigenstrat_20,
                                         n_iterations = "20")

iter_20_df<-iteration_20_results$pca_TAC_ld0.9_rm_out_liner %>% 
  dplyr::left_join(lineage_raw %>% select(isotype,Lineage), by = c("isotype"))

PC12_up<-iter_20_df %>% 
  dplyr::filter(PC2>0.04) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_left<-iter_20_df %>% 
  dplyr::filter(PC2<0.04) %>% 
  dplyr::filter(PC1<0.05) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_right<-iter_20_df %>% 
  dplyr::filter(PC1>0.2) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#### 25 iterations
TracyWidom_25<-"../../processed_data/PCA_more_iterations/PCA_25_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_25<-"../../processed_data/PCA_more_iterations/PCA_25_iterations/eigenstrat_outliers_removed.evac"

iteration_25_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_25,
                                         eigenstrat_file = eigenstrat_25,
                                         n_iterations = "25")

iter_25_df<-iteration_25_results$pca_TAC_ld0.9_rm_out_liner %>% 
  dplyr::left_join(lineage_raw %>% select(isotype,Lineage), by = c("isotype"))

PC12_up<-iter_25_df %>% 
  dplyr::filter(PC2>0.04) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_left<-iter_25_df %>% 
  dplyr::filter(PC2<0.04) %>% 
  dplyr::filter(PC1<0.05) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_right<-iter_25_df %>% 
  dplyr::ilter(PC1>0.2) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#### 30 iterations
TracyWidom_30<-"../../processed_data/PCA_more_iterations/PCA_30_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_30<-"../../processed_data/PCA_more_iterations/PCA_30_iterations/eigenstrat_outliers_removed.evac"

iteration_30_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_30,
                                         eigenstrat_file = eigenstrat_30,
                                         n_iterations = "30")

iter_30_df<-iteration_30_results$pca_TAC_ld0.9_rm_out_liner %>% 
  dplyr::left_join(lineage_raw %>% select(isotype,Lineage), by = c("isotype"))

PC12_up<-iter_30_df %>% 
  dplyr::filter(PC2>0.04) %>% 
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_left<-iter_30_df %>% 
  dplyr::filter(PC2<0.04) %>% 
  dplyr::filter(PC1<0.05) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

PC12_down_right<-iter_30_df %>% 
  dplyr::ilter(PC1>0.2) %>%
  dplyr::group_by(Lineage) %>% 
  dplyr::summarise(n=n())

#### 35 iterations
TracyWidom_35<-"../../processed_data/PCA_more_iterations/PCA_35_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_35<-"../../processed_data/PCA_more_iterations/PCA_35_iterations/eigenstrat_outliers_removed.evac"

iteration_35_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_35,
                                            eigenstrat_file = eigenstrat_35,
                                            n_iterations = "35")

#### 40 iterations
TracyWidom_40<-"../../processed_data/PCA_more_iterations/PCA_40_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_40<-"../../processed_data/PCA_more_iterations/PCA_40_iterations/eigenstrat_outliers_removed.evac"

iteration_40_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_40,
                                            eigenstrat_file = eigenstrat_40,
                                            n_iterations = "40")

######### assemble outliner removal plots #########
row1 <- cowplot::plot_grid(
  iteration_5_results$p_iteration, 
  iteration_10_results$p_iteration,
  ncol = 2,
  align = "h",
  axis = "tb"
)

row2 <- cowplot::plot_grid(
  iteration_15_results$p_iteration, 
  iteration_20_results$p_iteration,
  ncol = 2,
  align = "h",
  axis = "tb"
)

row3 <- cowplot::plot_grid(
  iteration_25_results$p_iteration, 
  iteration_30_results$p_iteration,
  ncol = 2,
  align = "h",
  axis = "tb"
)

final_plot <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  iteration_5_results$common_legend,
  ncol = 1,
  rel_heights = c(1, 1, 1, 
                  # 1,
                  0.2),  
  align = "v",
  axis = "lr"
)

ggsave(
  "../../figures/SF3_Cb_iteration_5_30_PCA.pdf",
  plot = final_plot,
  width = 7,
  height = 7,
  units = "in"
)

####################################
###### How many isotypes left ######
####################################

iteration_25_results$n_of_isotype
# 487
iteration_30_results$n_of_isotype
# 484
iteration_35_results$n_of_isotype
# 484
iteration_40_results$n_of_isotype
# 484
