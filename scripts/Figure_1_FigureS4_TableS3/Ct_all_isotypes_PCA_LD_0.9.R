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


geo_info_raw<-read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw 

# tracy_for_plot<-data.table::fread("../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") 


# calculate % var explained.
tracy_for_plot <- data.table::fread("../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") %>%
  dplyr::mutate(sum = sum(eigenvalue),
                VarExp = eigenvalue/sum,
                sigEV = ifelse(`p-value` < 0.05, T, F)) %>%
  dplyr::group_by(sigEV) %>%
  dplyr::mutate(sigVarExp = sum(VarExp)) %>%
  dplyr::ungroup()

# make priciple components df with annotated collections 0.9
cdf9no <- data.table::fread("../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/eigenstrat_no_removal.evac", skip = 1)

cdf9_2no <- cdf9no %>%
  dplyr::select(isotype = V1, PC1=V2, PC2=V3, PC3=V4, PC4=V5, PC5=V6, PC6 = V7) #  six significant eigenvectors wihtout outlier removal

# perform  clustering following https://www.datacamp.com/community/tutorials/hierarchical-clustering-R#what
# remove isotype to make distance matrix
cdf9_3no <- cdf9_2no %>% dplyr::select(-isotype)




# add back location data nd use these dfs to plot PC1 by PC2
cdf9_7no <- cdf9_2no %>%
  dplyr::left_join(geo_info,by=c("isotype")) 


# # export populations tmp 
# save(cdf9_7, cdf9_7no, file = "../../processed_data/Hawaii/cluster_assignments/tmp_cluster_assignments_LD0.9.Rdata")
# 

# Make named dataframe
pca_TAC_ld0.9_no_rm <- cdf9_7no 



###### Plot PCA function ######
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
    geom_point(shape=16, alpha=0.8, size=1.5, aes(x=.data[[x_axis]],y=.data[[y_axis]], color=geo))+
    scale_color_manual(values = geo.colours, name = "geo") +
    theme_bw() +
    labs(x=paste(x_axis," (",label_x,"%)",sep = ""),
         y=paste(y_axis," (",label_y,"%)",sep = ""))+
    theme(axis.title = element_text(size=8, face = "bold",color = "black"),
          axis.text = element_text(size=7, color = "black"),
          legend.position='none',
          # axis.title.x=element_blank(),
          panel.grid = element_blank())
  # labs(x=x_axis, y=y_axis)  
  # scale_x_continuous(breaks = c(0, 0.1, 0.2)) +
  # guides(fill= guide_legend(nrow=2))
  p1_1
  
  
  
  return(p1_1)
  
}


p_PC1_PC2<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC1",y_axis="PC2")
p_PC1_PC2
# ggsave("Ct_all_isotypes_PCA_plot_PC1_2_with_density.pdf", plot = p_PC1_PC2, width = 3.75, height = 3.75, units = "in")


p_PC3_PC4<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC3",y_axis="PC4")
p_PC3_PC4
# ggsave("Ct_all_isotypes_PCA_plot_PC3_4_with_density.pdf", plot = p_PC3_PC4, width = 3.75, height = 3.75, units = "in")





p_PC1_PC3<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC1",y_axis="PC3")
p_PC1_PC3



p_PC2_PC3<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC2",y_axis="PC3")
p_PC2_PC3






#### output PCA table
write.table(pca_TAC_ld0.9_no_rm,
            "Ct_pca_table_ld0.9.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)





library(patchwork)
Ct_all_isotypes_PCA_plot_PC1_4 <- (
  p_PC1_PC2 | p_PC3_PC4  
)
Ct_all_isotypes_PCA_plot_PC1_4

saveRDS(Ct_all_isotypes_PCA_plot_PC1_4, 
        file = "../../processed_data/assemble_figure_1/Ct_all_isotypes_PCA_plot_PC1_4.rds")

ggsave("Ct_all_isotypes_PCA_plot_PC1_4.pdf",
       Ct_all_isotypes_PCA_plot_PC1_4,
       height = 2.5, width = 7.5)




####################################
########## some details ###########
####################################
big_cluster<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC1 > -0.1) 
nrow(big_cluster)


left_up_cluster<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC1 < -0.1) %>% 
  filter(PC2 > 0) 
nrow(left_up_cluster)


left_down_cluster<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC1 < -0.1) %>% 
  filter(PC2 < 0) 
nrow(left_down_cluster)










