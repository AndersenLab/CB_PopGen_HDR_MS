rm(list=ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(dendextend)
library(ggpubr)
library(cowplot)
library(patchwork)
library(ggridges)


# source
source("../utilities.R")


geo_info_raw<-read.csv("../../processed_data/geo_info/Cb_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw 


# ###### Plot PCA function ######
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



########## remove out-liner iterations plot #######


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

# perform  clustering following https://www.datacamp.com/community/tutorials/hierarchical-clustering-R#what
# remove isotype to make distance matrix
cdf9_3no_rm_out_liner <- cdf9_2no_rm_out_liner %>% dplyr::select(-isotype)




# add back location data nd use these dfs to plot PC1 by PC2
cdf9_7no_rm_out_liner <- cdf9_2no_rm_out_liner %>%
  dplyr::left_join(geo_info,by=c("isotype")) 


# # export populations tmp 
# save(cdf9_7, cdf9_7no, file = "../../processed_data/Hawaii/cluster_assignments/tmp_cluster_assignments_LD0.9.Rdata")
# 

# Make named dataframe
pca_TAC_ld0.9_rm_out_liner <- cdf9_7no_rm_out_liner
n_of_isotype<-nrow(pca_TAC_ld0.9_rm_out_liner)
n_of_isotype
# 596



p_PC1_PC2_rm_out_liner<- plot_PCA(pca_TAC_ld0.9_rm_out_liner,
                                  tracy_for_plot_input = tracy_for_plot_rm_out_liner,
                                  x_axis="PC1",y_axis="PC2")
p_PC1_PC2_rm_out_liner
# ggsave("Cb_rm_out_liner_PCA_plot_PC1_2_with_density.pdf", plot = p_PC1_PC2_rm_out_liner, width = 3.75, height = 3.75, units = "in")


p_PC3_PC4_rm_out_liner<- plot_PCA(PCA_input=pca_TAC_ld0.9_rm_out_liner,
                                  tracy_for_plot_input = tracy_for_plot_rm_out_liner,
                                  x_axis="PC3",y_axis="PC4")
p_PC3_PC4_rm_out_liner
# ggsave("Cb_rm_out_liner_PCA_plot_PC3_4_with_density.pdf", plot = p_PC3_PC4, width = 3.75, height = 3.75, units = "in")



# generate legend
legend_plot_rm_out_liner <- ggplot(pca_TAC_ld0.9_rm_out_liner)+
  geom_point(aes(x = PC1, y = PC2, color = geo), alpha=0.8, size=1.2, shape = 16) +
  scale_color_manual(values = geo.colours, name = "Geo") +
  theme_bw() +
  # guides(color = guide_legend(nrow = 2, byrow = TRUE, keyheight = unit(0.01, "cm"))) +
  # theme(legend.position = "bottom",
  #       legend.spacing.y = unit(0, "cm"))
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

legend_plot_rm_out_liner


library(cowplot)
# legend_density <- get_legend(legend_plot)
legend_density_rm_out_liner <-cowplot::get_plot_component(legend_plot_rm_out_liner, "guide-box-bottom")
# 
# library(patchwork)
# Cb_rm_out_liner_PCA_plot_PC1_4_with_density <- (
#   p_PC1_PC2_rm_out_liner | p_PC3_PC4_rm_out_liner  
# ) / legend_density_rm_out_liner +              
#   plot_layout(heights = c(8.5, 1.5)) 



library(patchwork)
Cb_rm_out_liner_PCA_plot_PC1_4_with_density <- (
  p_PC1_PC2_rm_out_liner | p_PC3_PC4_rm_out_liner  
) 



Cb_rm_out_liner_PCA_plot_PC1_4_with_density<-Cb_rm_out_liner_PCA_plot_PC1_4_with_density +
  plot_annotation(
    title = paste0("PCA after ",
                   n_iterations,
                   " iterations of outliner removal  -  n_isotypes = ",
                   n_of_isotype)
  ) &
  theme(
    plot.title = element_text(size = 7, face = "bold", 
                              hjust = 0.5, vjust = -4),
  )

Cb_rm_out_liner_PCA_plot_PC1_4_with_density


return(list(p_iteration=Cb_rm_out_liner_PCA_plot_PC1_4_with_density,
            common_legend=legend_density_rm_out_liner,
            n_of_isotype=n_of_isotype
       )
)

}


#### 5 iterations
TracyWidom_5<-"../../processed_data/PCA_more_iterations/PCA_5_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_5<-"../../processed_data/PCA_more_iterations/PCA_5_iterations/eigenstrat_outliers_removed.evac"

iteration_5_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_5,
                                        eigenstrat_file = eigenstrat_5,
                                        n_iterations = "5")
iteration_5_results$p_iteration
ggsave(
  "Cb_rm_out_liner_5_iterations_PCA_plot_PC1_4_with_density.pdf", 
  plot = iteration_5_results$p_iteration, 
  width = 7, 
  height = 3,   
  units = "in"
)







#### 10 iterations
TracyWidom_10<-"../../processed_data/PCA_more_iterations/PCA_10_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_10<-"../../processed_data/PCA_more_iterations/PCA_10_iterations/eigenstrat_outliers_removed.evac"

iteration_10_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_10,
                                        eigenstrat_file = eigenstrat_10,
                                        n_iterations = "10")
iteration_10_results$p_iteration






#### 15 iterations
TracyWidom_15<-"../../processed_data/PCA_more_iterations/PCA_15_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_15<-"../../processed_data/PCA_more_iterations/PCA_15_iterations/eigenstrat_outliers_removed.evac"

iteration_15_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_15,
                                         eigenstrat_file = eigenstrat_15,
                                         n_iterations = "15")
iteration_15_results$p_iteration




#### 20 iterations
TracyWidom_20<-"../../processed_data/PCA_more_iterations/PCA_20_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_20<-"../../processed_data/PCA_more_iterations/PCA_20_iterations/eigenstrat_outliers_removed.evac"

iteration_20_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_20,
                                         eigenstrat_file = eigenstrat_20,
                                         n_iterations = "20")
iteration_20_results$p_iteration







#### 25 iterations
TracyWidom_25<-"../../processed_data/PCA_more_iterations/PCA_25_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_25<-"../../processed_data/PCA_more_iterations/PCA_25_iterations/eigenstrat_outliers_removed.evac"

iteration_25_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_25,
                                         eigenstrat_file = eigenstrat_25,
                                         n_iterations = "25")
iteration_25_results$p_iteration





#### 30 iterations
TracyWidom_30<-"../../processed_data/PCA_more_iterations/PCA_30_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_30<-"../../processed_data/PCA_more_iterations/PCA_30_iterations/eigenstrat_outliers_removed.evac"

iteration_30_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_30,
                                         eigenstrat_file = eigenstrat_30,
                                         n_iterations = "30")
iteration_30_results$p_iteration







#### 35 iterations
TracyWidom_35<-"../../processed_data/PCA_more_iterations/PCA_35_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_35<-"../../processed_data/PCA_more_iterations/PCA_35_iterations/eigenstrat_outliers_removed.evac"

iteration_35_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_35,
                                            eigenstrat_file = eigenstrat_35,
                                            n_iterations = "35")
iteration_35_results$p_iteration










#### 40 iterations
TracyWidom_40<-"../../processed_data/PCA_more_iterations/PCA_40_iterations/TracyWidom_statistics_outlier_removal.tsv"
eigenstrat_40<-"../../processed_data/PCA_more_iterations/PCA_40_iterations/eigenstrat_outliers_removed.evac"

iteration_40_results<-plot_remove_out_liner(TracyWidom_file=TracyWidom_40,
                                            eigenstrat_file = eigenstrat_40,
                                            n_iterations = "40")
iteration_40_results$p_iteration







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


row4 <- cowplot::plot_grid(
  iteration_35_results$p_iteration, 
  iteration_40_results$p_iteration,
  ncol = 2,
  align = "h",
  axis = "tb"
)



final_plot <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  row4,
  iteration_5_results$common_legend,
  ncol = 1,
  rel_heights = c(1, 1, 1, 1, 0.2),  
  align = "v",
  axis = "lr"
)


ggsave(
  "Cb_iteration_5_30_PCA_cowplot.pdf",
  plot = final_plot,
  width = 7,
  height = 7,
  units = "in"
)






####################################
###### How many isotypes left ######
####################################

iteration_25_results$n_of_isotype
# 410
iteration_30_results$n_of_isotype
# 408
iteration_35_results$n_of_isotype
# 408
iteration_40_results$n_of_isotype
# 408




