rm(list=ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(GGally)
library(dendextend)
library(ggpubr)
library(cowplot)
library(patchwork)
library(ggridges)

source("../utilities.R")

geo_info_raw <- read.csv("../../processed_data/geo_info/Cb_indep_isotype_info_geo.csv")
geo_info <- geo_info_raw

plot_PCA <- function(PCA_input, tracy_for_plot_input, x_axis, y_axis){
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
    theme(axis.title = element_text(size=7, face = "bold",color = "black"),
          axis.text = element_text(size=6, color = "black"),
          legend.position='none',
          panel.grid = element_blank())
  
  return(p1_1)
}

base_dir <- "../../processed_data/PCA_by_chrom"

chroms <- c("I","II","III","IV","V","X")

combined_panels <- list()
pcs_geo_by_chrom <- list()

for (chrom in chroms) {
  message("Processing chromosome: ", chrom)
  
  tracy_path <- file.path(base_dir, chrom, "EIGESTRAT", "LD_0.9", "NO_REMOVAL", "TracyWidom_statistics_no_removal.tsv")
  evac_path  <- file.path(base_dir, chrom, "EIGESTRAT", "LD_0.9", "NO_REMOVAL", "eigenstrat_no_removal.evac")
  
  if (!file.exists(tracy_path)) {
    warning("Tracy file not found for ", chrom, " at: ", tracy_path, "  -> skipping this chromosome.")
    next
  }
  if (!file.exists(evac_path)) {
    warning("Evac file not found for ", chrom, " at: ", evac_path, "  -> skipping this chromosome.")
    next
  }
  
  tracy_for_plot <- data.table::fread(tracy_path) %>%
    dplyr::mutate(sum = sum(eigenvalue),
                  VarExp = eigenvalue / sum,
                  sigEV = ifelse(`p-value` < 0.05, TRUE, FALSE)) %>%
    dplyr::group_by(sigEV) %>%
    dplyr::mutate(sigVarExp = sum(VarExp)) %>%
    dplyr::ungroup()
  
  
  
  
  
  cdf <- data.table::fread(evac_path, skip = 1)
  
  cdf_pcs <- cdf %>%
    dplyr::select(isotype = V1, PC1 = V2, PC2 = V3, PC3 = V4, PC4 = V5, PC5 = V6, PC6 = V7)
  
  cdf_pcs_geo <- cdf_pcs %>%
    dplyr::left_join(geo_info, by = c("isotype"))
  
  pcs_geo_by_chrom[[chrom]] <- cdf_pcs_geo
  
  
  
  
  
  p_PC1_PC2 <- plot_PCA(PCA_input = cdf_pcs_geo,
                        tracy_for_plot_input = tracy_for_plot,
                        x_axis = "PC1", y_axis = "PC2")
  
  p_PC3_PC4 <- plot_PCA(PCA_input = cdf_pcs_geo,
                        tracy_for_plot_input = tracy_for_plot,
                        x_axis = "PC3", y_axis = "PC4")
  
  
  combined_core <- cowplot::plot_grid(p_PC1_PC2, p_PC3_PC4, ncol = 2, align = "hv", rel_heights = c(1,1))
  
  combined_titled <- cowplot::ggdraw() +
    cowplot::draw_plot(combined_core) +
    cowplot::draw_label(paste0("Chromosome ", chrom),
                        x = 0.5, y = 0.97, hjust = 0.5, vjust = 0,
                        fontface = "bold", size = 8)
  
  combined_panels[[chrom]] <- combined_titled
  
}

if (length(combined_panels) == 0) stop("No panels were created (missing files or all skipped).")

final_plot <- patchwork::wrap_plots(plotlist = combined_panels, ncol = 2)


ggsave("PCA_by_chrom_raw.pdf", final_plot, width = 7, height = 4, units = "in", device = "pdf")





library(readr)
lineage_raw<-readr::read_tsv("../../data/From_Nic/before_20250828_isotype_byLineage_GeoLocAdmCol.tsv")
lineage<-lineage_raw %>% 
  select(isotype,Lineage)


chr_num<-c(1:6)
for (i in chr_num){
  print(
    pcs_geo_by_chrom[[i]] %>% 
      left_join(lineage,by = c("isotype")) %>% 
      filter(PC1 >0.1) %>% 
      count(Lineage)
  )
}




ggplot(pcs_geo_by_chrom[[1]] %>% filter(PC1 <0.05), aes(x = PC2)) +
  geom_histogram(bins = 100)

View(pcs_geo_by_chrom[[1]])

ChrI_PC12_mixed<-pcs_geo_by_chrom[[1]]%>% 
  filter(PC1 <0.05) %>% 
  filter(PC2 <0) %>% 
  left_join(lineage, by = c("isotype"))

write.csv(ChrI_PC12_mixed,
          "ChrI_PC12_mixed.csv",
          row.names = FALSE,
          quote = FALSE)



pcs_geo_by_chrom[[1]]%>% 
  filter(PC4 < -0.2) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[1]]%>% 
  filter(PC3 > 0.1) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[1]]%>% 
  filter(PC4 > 0.05) %>% 
  left_join(lineage, by = c("isotype"))



ggplot(pcs_geo_by_chrom[[1]] %>% filter(PC4 > -0.1 & PC4 < 0.05), aes(x = PC3)) +
  geom_histogram(bins = 100)

ChrI_PC34_mixed<-pcs_geo_by_chrom[[1]] %>% 
  filter(PC4 > -0.1 & PC4 < 0.05) %>% 
  filter(PC3 < 0 & PC3 > -0.075)












ggplot(pcs_geo_by_chrom[[2]] %>% filter(PC1 <0.05), aes(x = PC2)) +
  geom_histogram(bins = 100)

View(pcs_geo_by_chrom[[2]])

ChrII_PC12_mixed<-pcs_geo_by_chrom[[2]]%>% 
  filter(PC1 <0.05) %>% 
  filter(PC2 <0) %>% 
  left_join(lineage, by = c("isotype"))

write.csv(ChrII_PC12_mixed,
          "ChrII_PC12_mixed.csv",
          row.names = FALSE,
          quote = FALSE)



pcs_geo_by_chrom[[2]]%>% 
  filter(PC4 > 0.4) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[2]]%>% 
  filter(PC3 > 0.2) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[2]]%>% 
  filter(PC3 < -0.03) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[2]]%>% 
  filter(PC4 < -0.05) %>% 
  left_join(lineage, by = c("isotype"))



ggplot(pcs_geo_by_chrom[[2]] %>% filter(PC4 > -0.05 & PC4 < 0.2 & PC3 < 0.2 & PC3 > -0.03), aes(x = PC4)) +
  geom_histogram(bins = 100)

ChrII_PC34_mixed<-pcs_geo_by_chrom[[2]] %>%
  filter(PC4 > 0.05 & PC4 < 0.2)

ChrII_PC34_global<-pcs_geo_by_chrom[[2]] %>%
  filter(PC4 > -0.05 & PC4 < 0.2 & PC3 < 0.2 & PC3 > -0.03) %>% 
  filter(PC4 < 0.05)













ggplot(pcs_geo_by_chrom[[3]] %>% filter(PC1 <0.05), aes(x = PC2)) +
  geom_histogram(bins = 100)

View(pcs_geo_by_chrom[[3]])

ChrIII_PC12_global<-pcs_geo_by_chrom[[3]]%>% 
  filter(PC1 < 0.05) %>% 
  filter(PC2 > -0.05) %>% 
  left_join(lineage, by = c("isotype"))

ChrIII_PC12_mixed<-pcs_geo_by_chrom[[3]]%>% 
  filter(PC1 < 0.05) %>% 
  filter(PC2 < -0.05) %>% 
  filter(PC2 > -0.15) %>% 
  left_join(lineage, by = c("isotype"))

write.csv(ChrIII_PC12_mixed,
          "ChrIII_PC12_mixed.csv",
          row.names = FALSE,
          quote = FALSE)



pcs_geo_by_chrom[[3]]%>% 
  filter(PC3 > 0.1) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[3]]%>% 
  filter(PC4 < -0.3) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[3]]%>% 
  filter(PC3 < 0.1) %>% 
  filter(PC4 > -0.3) %>% 
  filter(PC4 < -0.05) %>% 
  left_join(lineage, by = c("isotype"))

pcs_geo_by_chrom[[3]]%>% 
  filter(PC3 < 0.1) %>% 
  filter(PC4 > -0.3) %>% 
  filter(PC4 > -0.05) %>% 
  left_join(lineage, by = c("isotype"))








ggplot(pcs_geo_by_chrom[[4]] %>% filter(PC4 > -0.05 & PC4 < 0.05), aes(x = PC3)) +
  geom_histogram(bins = 100)

ggplot(pcs_geo_by_chrom[[4]], aes(x = PC4)) +
  geom_histogram(bins = 100)


View(pcs_geo_by_chrom[[4]])

ChrIV_PC34_global<-pcs_geo_by_chrom[[4]]%>% 
  filter(PC3 > 0) %>% 
  filter(PC3 < 0.025) %>% 
  left_join(lineage, by = c("isotype"))

ChrIV_PC34_mixed_3<-pcs_geo_by_chrom[[4]]%>% 
  filter(PC3 > -0.05) %>% 
  filter(PC3 < 0) %>% 
  left_join(lineage, by = c("isotype"))


ChrIV_PC34_mixed_1<-pcs_geo_by_chrom[[4]]%>% 
  filter(PC4 > 0.05) %>% 
  left_join(lineage, by = c("isotype"))


ChrIV_PC34_mixed_2<-pcs_geo_by_chrom[[4]]%>% 
  filter(PC4 > -0.05 & PC4 < 0.05) %>% 
  filter(PC3 < -0.05) %>% 
  left_join(lineage, by = c("isotype"))










ggplot(pcs_geo_by_chrom[[5]] %>% filter(PC4 > -0.05 & PC4 < 0.05), aes(x = PC3)) +
  geom_histogram(bins = 100)


View(pcs_geo_by_chrom[[5]])

ChrV_PC34_Taiwan<-pcs_geo_by_chrom[[5]]%>% 
  filter(PC4 < -0.1) %>% 
  left_join(lineage, by = c("isotype"))

ChrV_PC34_global<-pcs_geo_by_chrom[[5]]%>% 
  filter(PC3 > -0.05) %>% 
  left_join(lineage, by = c("isotype"))

ChrV_PC34_global<-pcs_geo_by_chrom[[5]]%>% 
  filter(PC4 < 0.12 & PC4 > 0.01) %>% 
  left_join(lineage, by = c("isotype"))










ggplot(pcs_geo_by_chrom[[6]] %>% filter(PC1 < 0.05), aes(x = PC2)) +
  geom_histogram(bins = 100)


View(pcs_geo_by_chrom[[6]])

ChrX_PC12_global<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC1 < 0.05) %>% 
  filter(PC2 > -0.025) %>% 
  left_join(lineage, by = c("isotype"))

ChrX_PC12_mixed<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC1 < 0.05) %>% 
  filter(PC2 < -0.025) %>% 
  left_join(lineage, by = c("isotype"))



ChrX_PC34_India<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC3 > 0.2) %>% 
  left_join(lineage, by = c("isotype"))



ChrX_PC34_Indo_and_CosmopilitanNairobi<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC3 < 0.2 & PC4 > 0.27) %>% 
  left_join(lineage, by = c("isotype"))

ChrX_PC34_middle<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC3 < 0.2 & PC4 > 0.15 & PC4 < 0.27) %>% 
  left_join(lineage, by = c("isotype"))


ChrX_PC34_Taiwan<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC3 < 0.2 & PC4 < -0.05) %>% 
  left_join(lineage, by = c("isotype"))



ChrX_PC34_Australia<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC3 < 0.2 & PC4 > -0.05 & PC4 < 0.15 & PC3 < -0.01) %>% 
  left_join(lineage, by = c("isotype"))


ChrX_PC34_Australia<-pcs_geo_by_chrom[[6]]%>% 
  filter(PC3 < 0.2 & PC4 > -0.05 & PC4 < 0.15 & PC3 > -0.01) %>% 
  left_join(lineage, by = c("isotype"))

