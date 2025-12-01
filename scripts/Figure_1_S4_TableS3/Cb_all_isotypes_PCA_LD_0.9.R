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

# tracy_for_plot<-data.table::fread("../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") 


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
    theme(axis.title = element_text(size=10, face = "bold",color = "black"),
          axis.text = element_text(size=9, color = "black"),
          legend.position='none',
          # legend.position='bottom',
          # axis.title.x=element_blank(),
          panel.grid = element_blank())
  
  return(p1_1)
 
  
}


p_PC1_PC2<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC1",y_axis="PC2")
p_PC1_PC2



p_PC3_PC4<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC3",y_axis="PC4")
p_PC3_PC4




# 
# #### output PCA table
# write.table(pca_TAC_ld0.9_no_rm,
#             "Cb_pca_table_ld0.9.tsv",
#             sep = '\t',
#             col.names = TRUE,
#             row.names = FALSE,
#             quote = FALSE)
# 




library(patchwork)
Cb_all_isotypes_PCA_plot_PC1_4 <- (
  p_PC1_PC2 | p_PC3_PC4
) + plot_layout(guides = "collect") & 
  scale_color_manual(values = geo.colours, name = "geo",breaks = "Cosmopolitan") &
  # scale_color_discrete(breaks = "Cosmopolitan") & 
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


Cb_all_isotypes_PCA_plot_PC1_4

saveRDS(Cb_all_isotypes_PCA_plot_PC1_4, 
        file = "../../processed_data/assemble_figure_1/Cb_all_isotypes_PCA_plot_PC1_4.rds")

saveRDS(p_PC1_PC2$legend, 
        file = "../../processed_data/assemble_figure_1/Cb_PCA_legend.rds")



# ggsave("Cb_all_isotypes_PCA_plot_PC1_4.pdf",
#        Cb_all_isotypes_PCA_plot_PC1_4,
#        height = 3, width = 5.25)







library(dplyr)
PC12KD<-pca_TAC_ld0.9_no_rm %>%
  filter(PC1 > 0.1 & PC2 < -0.2)
nrow(PC12KD)
# 10 KD


PC12AD<-pca_TAC_ld0.9_no_rm %>%
  filter(PC1 > 0.1 & PC2 > 0)
nrow(PC12AD)
#32 AD


ggplot(pca_TAC_ld0.9_no_rm, aes(x = PC3)) +
  geom_histogram()

ggplot(pca_TAC_ld0.9_no_rm, aes(x = PC4)) +
  geom_histogram(bins = 100)


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

# write.csv(PC34mix,
#           "PC34mix.csv",
#           row.names = FALSE,
#           quote = FALSE)
#45
### add previous lineage 
library(readr)
library(dplyr)
lineage_raw<-read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byLineage_GeoLocAdmCol_20250909.tsv")
# lineage_raw<-read_tsv("../../data/From_Nic/isotype_byLineage_GeoLocAdmCol_20250828.tsv")
lineage<-lineage_raw %>% select(isotype,Lineage)
PC34mix<-PC34mix %>% 
  left_join(lineage)

# write.csv(PC34mix,
#           "PC34mix.csv",
#           row.names = FALSE,
#           quote = FALSE)






PC34Indo<-pca_TAC_ld0.9_no_rm %>%
  filter(PC4 < -0.4)
nrow(PC34Indo)






########### info #######
tracy_for_plot2 <- tracy_for_plot %>%
  arrange(N) %>% 
  mutate(cumVarExp = cumsum(VarExp)) 









####### Table S3 ####

### PCA differentiated isotype clusters. 

PCA_differentiated<-rbind((PC12AD %>% mutate(Cluster = c("Australia")) ),
                          
                          (PC12KD%>% mutate(Cluster = c("India")) ),
                          
                          (PC34mix %>% select(-Lineage) %>% mutate(Cluster = c("Mixed")) ),
                          
                          (PC34TW %>% mutate(Cluster = c("Taiwan")) ),
                          
                          (PC34Indo %>% mutate(Cluster = c("Indonesia")) )
                          
) 

#### output PCA table S3 table
TableS3<-pca_TAC_ld0.9_no_rm %>% 
  select(-PC5,-PC6,-lat,-long) %>% 
  left_join(PCA_differentiated %>% select(isotype,Cluster)) %>% 
  mutate(Cluster = ifelse(is.na(Cluster),"Global",Cluster))

write.table(TableS3,
            "../../tables/TableS3_PCA_data.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)





    