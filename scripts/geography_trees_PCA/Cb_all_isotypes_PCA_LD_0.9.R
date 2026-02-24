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

# calculate % var explained.
tracy_for_plot <- data.table::fread("../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") %>%
  dplyr::mutate(sum = sum(eigenvalue),
                VarExp = eigenvalue/sum,
                sigEV = ifelse(`p-value` < 0.05, T, F)) %>%
  dplyr::group_by(sigEV) %>%
  dplyr::mutate(sigVarExp = sum(VarExp)) %>%
  dplyr::ungroup()

# make priciple components df with annotated collections 0.9
cdf9no <- data.table::fread("../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/eigenstrat_no_removal.evac", skip = 1)

cdf9_2no <- cdf9no %>%
  dplyr::select(isotype = V1, PC1=V2, PC2=V3, PC3=V4, PC4=V5, PC5=V6, PC6 = V7) #  six significant eigenvectors wihtout outlier removal

cdf9_3no <- cdf9_2no %>% dplyr::select(-isotype)


# add back location data nd use these dfs to plot PC1 by PC2
cdf9_7no <- cdf9_2no %>%
  dplyr::left_join(geo_info,by=c("isotype")) 

# Make named dataframe
pca_TAC_ld0.9_no_rm <- cdf9_7no 


###### Plot PCA function ######
plot_PCA<-function(PCA_input, tracy_for_plot_input, x_axis, y_axis){
  
  label_x<-tracy_for_plot_input %>% 
    dplyr::filter(N == as.numeric(sub("PC", "", x_axis))) %>% 
    dplyr::select(VarExp) %>%
    dplyr::pull(VarExp) %>%
    `*`(100) %>%
    round(digits = 2) %>%
    as.character()
  
  label_y<-tracy_for_plot_input %>% 
    dplyr::filter(N == as.numeric(sub("PC", "", y_axis))) %>% 
    dplyr::select(VarExp) %>%
    dplyr::pull(VarExp) %>%
    `*`(100) %>%
    round(digits = 2) %>%
    as.character()
  
  p1_1<-ggplot2::ggplot(PCA_input)+
    geom_point(shape=16, alpha=0.8, size=1.5, aes(x=.data[[x_axis]],y=.data[[y_axis]], color=geo))+
    scale_color_manual(values = geo.colours, name = "geo") +
    theme_bw() +
    labs(x=paste(x_axis," (",label_x,"%)",sep = ""),
         y=paste(y_axis," (",label_y,"%)",sep = ""))+
    theme(axis.title = element_text(size=10, face = "bold",color = "black"),
          axis.text = element_text(size=9, color = "black"),
          legend.position='none',
          panel.grid = element_blank())
  
  return(p1_1)
 
}


p_PC1_PC2<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC1",y_axis="PC2")

p_PC3_PC4<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC3",y_axis="PC4")

Cb_all_isotypes_PCA_plot_PC1_4 <- (
  p_PC1_PC2 | p_PC3_PC4
) + plot_layout(guides = "collect") & 
  scale_color_manual(values = geo.colours, name = "geo",breaks = "Cosmopolitan") &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(0.3, "lines"),
    legend.text = element_text(size = 9),
    legend.spacing.x = unit(0.1, "cm"),
    legend.spacing.y = unit(0.1, "cm")
  ) & 
  guides(color = guide_legend(title = NULL, 
                              ncol = 6, 
                              nrow = 3,
                              byrow = TRUE))

saveRDS(Cb_all_isotypes_PCA_plot_PC1_4, 
        file = "../../processed_data/Geo_info/Cb_all_isotypes_PCA_plot_PC1_4.rds")

PC12KD<-pca_TAC_ld0.9_no_rm %>%
  filter(PC1 > 0.1 & PC2 < -0.2)
nrow(PC12KD)
# 10 KD

PC12AD<-pca_TAC_ld0.9_no_rm %>%
  filter(PC1 > 0.1 & PC2 > 0)
nrow(PC12AD)
#32 AD

PC34KD<-pca_TAC_ld0.9_no_rm %>%
  filter(PC3 > 0.07)
nrow(PC34KD)
#10

PC34TW<-pca_TAC_ld0.9_no_rm %>%
  filter(PC3 < -0.05 & PC3 > -0.15) %>% 
  filter(PC4 >0.06)
nrow(PC34TW)
#27

PC34mix<-pca_TAC_ld0.9_no_rm %>%
  filter(PC3 < -0.05) %>% 
  filter(PC4 <0.06 & PC4 > -0.2)
nrow(PC34mix)

### add previous lineage 
lineage_raw<-read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv")
lineage<-lineage_raw %>% select(isotype,Lineage)
PC34mix<-PC34mix %>% 
  dplyr::left_join(lineage)

PC34Indo<-pca_TAC_ld0.9_no_rm %>%
  dplyr::filter(PC4 < -0.4)
nrow(PC34Indo)


########### info #######
tracy_for_plot2 <- tracy_for_plot %>%
  dplyr::arrange(N) %>% 
  dplyr::mutate(cumVarExp = cumsum(VarExp)) 


####### Table S3 ####
### PCA differentiated isotype clusters. 
PCA_differentiated<-rbind((PC12AD %>% mutate(Cluster = c("Australia")) ),
                          (PC12KD%>% mutate(Cluster = c("India")) ),
                          (PC34mix %>% select(-Lineage) %>% mutate(Cluster = c("Mixed")) ),
                          (PC34TW %>% mutate(Cluster = c("Taiwan")) ),
                          (PC34Indo %>% mutate(Cluster = c("Indonesia")) )
) 

#### output PCA table S3 table
SD2<-pca_TAC_ld0.9_no_rm %>% 
  dplyr::select(-PC5,-PC6,-lat,-long) %>% 
  dplyr::left_join(PCA_differentiated %>% select(isotype,Cluster)) %>% 
  dplyr::mutate(Cluster = ifelse(is.na(Cluster),"Global",Cluster))

write.table(SD2,
            "../../supplementary_data/SD2_PCA_data.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

    