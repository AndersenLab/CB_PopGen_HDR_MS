
rm(list=ls())


library(pophelper)
library(tidyverse)
library(ggthemes)

source("../utilities.R")


# My colors 
ancestry.colours <- setNames(
  c("#4b0200","#da000f", "#ff6d93", "#d45700",
    "#563900", "#ffe5ca", "#ffb914", "#ffda90", "#77c000",
    "#01e51b", "#00491e", "#01dea2", "#82fffa", "#00a4b1",
    "#5cb9ff", "#000f2d", "#0141b9", "#9c87ff", "#cbb1ff",
    "#f479ff", "#5f0058", "#BAB465", "#ad0041","#ffb4a8","black","grey"
    
  ),
  c(LETTERS[1:26])
)





isotype_geo_info<-read.csv("../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv") 

sample_list<-read.table("../../processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt")
library(dplyr)

desired <- as.character(sample_list$V1)

isotype_geo_info_ordered <- isotype_geo_info %>%
  filter(isotype %in% desired) %>%
  mutate(isotype = factor(isotype, levels = desired)) %>%
  arrange(isotype)


# get list of isotype names 
sample_names <- isotype_geo_info_ordered$isotype




# plot panel A, K by CV summary
k_summary <- data.table::fread("../../processed_data/Cb_admixture/admix_replicates_CV.tsv",header = F) 
colnames(k_summary)<-c("K","CV")
numeric_levels <- sort(as.numeric(gsub("K=", "", k_summary$K)))
k_summary$K <- factor(k_summary$K, levels = unique(paste0("K=", numeric_levels)))

set.seed(123)
ksum_plot <- ggplot(k_summary)+
  aes(x = factor(K), y = CV)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1,size=0.5)+
  theme_bw()+
  labs(x = "K")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ksum_plot
ggsave(ksum_plot, filename = "Cb_ksum_plot.pdf", height = 7, width = 7)


best_k <- data.frame(K.22 = 20:24)
best_k

k_values <- best_k$K.22
best_k_value <- as.numeric(sub("K.", "", colnames(best_k)[1]))


which_replicate<-c(1:10)

for (which_replicate in which_replicate) {
  admix_plots <- list()
  
for(kpops in 1:length(grep(".Q", list.files("../../processed_data/Cb_admixture/"), value = T))){
  K <- as.numeric(sub(".Q", "", strsplit(grep(".Q", list.files("../../processed_data/Cb_admixture/"), value = T)[kpops], split = "\\_")[[1]][3]))
  
  if (!(K %in% k_values)) {
    next  
  }
  
  # load Q files
  qfile_name <- grep(pattern = glue::glue("_{K}_\\d+\\.Q$"), value = T, x = list.files("../../processed_data/Cb_admixture/"))
  qfile <- pophelper::readQ(files = paste0("../../processed_data/Cb_admixture/",qfile_name))[[which_replicate]]
  # add pop names
  names_pool <- c(LETTERS, paste0("A", LETTERS))
  colnames(qfile) <- names_pool[1:K]
  
  
  qfile <- qfile %>%
    dplyr::mutate(samples = sample_names)
  
  
  
  if (K == best_k_value) {
    write.table(qfile, file = paste0("../../processed_data/Cb_admixture/K", K, "_Processed_Ancestry_replicate_",which_replicate,".tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  
  long_admix_pops <- qfile %>%
    dplyr::mutate(samples = sample_names) %>%
    tidyr::gather(cluster, frac_cluster, -samples) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(max_frac = max(frac_cluster)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  plot_order <- long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  

    admix_plots[[kpops]] <-long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
    theme_bw() +
    theme(axis.text.x=element_blank(),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(size = 8, angle = 90, vjust = .5),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))

  if(!exists("representative_K_strains")){
    representative_K_strains <- dplyr::filter(plot_order, max_frac  > 0.999
                                              # , !samples %in% names(strain_islands)
                                              ) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(sample_n = 1:n()) %>%
      dplyr::top_n(3, sample_n) %>%
      dplyr::mutate(K_size = K)
  } else {
    representative_K_strains <- dplyr::filter(plot_order, frac_cluster > 0.999
                                              ) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(sample_n = 1:n()) %>%
      dplyr::top_n(3, sample_n) %>%
      dplyr::mutate(K_size = K) %>%
      dplyr::bind_rows(representative_K_strains, .)
  }

  admix_plots <- Filter(Negate(is.null), admix_plots)
  
}




samples <- isotype_geo_info_ordered[, "isotype"]
groups <- isotype_geo_info_ordered[, "geo"]

sample_colors <- geo.colours[groups]



admix_plots[[1]] <- admix_plots[[1]] + theme(legend.position = "none")
admix_plots[[11]] <- admix_plots[[11]] + theme(legend.position = "none")
admix_plots[[21]] <- admix_plots[[21]] + theme(legend.position = "none")
admix_plots[[31]] <- admix_plots[[31]] + theme(legend.position = "none")

admix_plots_legend <- admix_plots[[41]] +
guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom") 
legend <- cowplot::get_legend(admix_plots_legend)


admix_plots[[41]] <- admix_plots[[41]] + theme(legend.position = "none")  


admixture_plots <- cowplot::plot_grid(admix_plots[[1]],
                   admix_plots[[11]],
                   admix_plots[[21]],
                   admix_plots[[31]],
                   admix_plots[[41]],
                   ncol = 1)
admixture_plots

final_admixture_plots <- cowplot::plot_grid(
  admixture_plots,
  legend,
  ncol = 1,
  rel_heights = c(4, 0.3)
)
final_admixture_plots

ggsave(final_admixture_plots, 
       filename = paste0("Cb_5_admixture_replicate_",which_replicate,".pdf"), height = 7, width = 7, useDingbats=FALSE)




ksummary_plot <- cowplot::plot_grid(ksum_plot,
                                    final_admixture_plots,
                                    ncol = 1,
                                    nrow = 2,
                                    labels = c("a", "b"),
                                    rel_heights = c(2, 3),
                                    label_y = c(1, 1.05))
ksummary_plot

}


which_replicate<-c(1:10)

for (which_replicate in which_replicate){

best_k_qfile<-read.table(paste0("../../processed_data/Cb_admixture/K22_Processed_Ancestry_replicate_",which_replicate,".tsv"),
                          header = TRUE)

best_k_long_admix_pops <- best_k_qfile %>%
  dplyr::mutate(samples = sample_names) %>%
  tidyr::gather(cluster, frac_cluster, -samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>%
  left_join(isotype_geo_info_ordered, by = c("samples" = "isotype"))

best_k_plot_order <- best_k_long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))


write.csv(best_k_long_admix_pops,
          paste0("../../processed_data/Cb_admixture/K22_best_k_long_admix_pops_replicate_",which_replicate,".csv")
          )


best_k_admix_plots<-best_k_long_admix_pops %>%
  dplyr::mutate(ordered_samples = factor(samples, levels = best_k_plot_order$samples)) %>%
  ggplot() +
  geom_bar(stat = "identity", 
           aes(x = ordered_samples, 
               y = frac_cluster, 
               fill = cluster)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(fill = "", x = "", y =  glue::glue("K = {best_k_value}")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1,
                                   size = 1,
                                   margin = margin(t = 0)),    
        axis.text.y=element_blank(),
        axis.title.y = element_text(size = 8,angle = 90, vjust = .5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.02, "cm"),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  coord_cartesian(expand = FALSE) + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))+
  guides(fill = guide_legend(nrow = 2)) + 
  facet_grid(.~ geo, scales = "free_x", space = "free_x")+
  theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                    margin = margin(b = 1, t = 3)),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")
        )
  


best_k_admix_plots

if (which_replicate == 1) {
ggsave(best_k_admix_plots, filename = "best_k_admix_plots_non_cosmopolitan_rep_1.pdf", height = 2, width = 9)
}

pie_best_k_long_admix_pops <- best_k_long_admix_pops %>%
  group_by(geo, cluster) %>%
  summarise(total_frac = sum(frac_cluster), .groups = 'drop') %>%
  group_by(geo) %>%
  mutate(percent = total_frac / sum(total_frac) * 100) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))))

stacked_bar <- ggplot(pie_best_k_long_admix_pops, 
                      aes(x = geo, y = percent, fill = cluster, order = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = ancestry.colours) +
  labs(
    x = NULL,
    y = "Percentage (%)",
    fill = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  labs(title = "All isotypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 10))


stacked_bar













library(ggpubr)

combined_plot <- ggpubr::ggarrange(best_k_admix_plots, stacked_bar,
                           ncol = 1, nrow = 2, 
                           heights = c(1, 1))
combined_plot














best_k_long_admix_pops <- best_k_qfile %>%
  dplyr::mutate(samples = sample_names) %>%
  tidyr::gather(cluster, frac_cluster, -samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>%
  left_join(isotype_geo_info_ordered, by = c("samples" = "isotype"))

rep_strains_best_k_long_admix_pops<-best_k_long_admix_pops %>% 
  filter(max_frac > 0.999)

rep_strains_best_k_plot_order <- rep_strains_best_k_long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))



# plot
rep_strains_best_k_admix_plots<-rep_strains_best_k_long_admix_pops %>%
  dplyr::mutate(ordered_samples = factor(samples, levels = best_k_plot_order$samples)) %>%
  ggplot() +
  geom_bar(stat = "identity", 
           aes(x = ordered_samples, 
               y = frac_cluster, 
               fill = cluster)) +
  scale_fill_manual(values = ancestry.colours) +
  labs(fill = "", x = "", y =  glue::glue("K = {best_k_value}")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1,
                                   size = 1,
                                   margin = margin(t = 0)),    
        axis.text.y=element_blank(),
        axis.title.y = element_text(size = 8, angle = 90, vjust = .5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.02, "cm"),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  coord_cartesian(expand = FALSE) +  
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))+
  guides(fill = guide_legend(nrow = 2)) + 
  facet_grid(.~ geo, scales = "free_x", space = "free_x")+
  theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                    margin = margin(b = 1, t = 3)),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")
  )



rep_strains_best_k_admix_plots


non_admixed_isotype_list<-rep_strains_best_k_plot_order %>% 
  select(-frac_cluster,-max_frac)

write.table(non_admixed_isotype_list,
            paste0("../../processed_data/Geo_info/non_admixed_isotype_replicate_",which_replicate,".txt"),
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)


n_each_geo<-isotype_geo_info_ordered %>% 
  group_by(geo) %>%
  summarise(total = n())

admixed_in_each_geo<-isotype_geo_info_ordered %>% 
  filter(!isotype %in% non_admixed_isotype_list$samples) %>% 
  group_by(geo) %>%
  summarise(admixed = n()
    ) %>% 
  left_join(n_each_geo, by = "geo") %>% 
  mutate(admixed_percentage = admixed / total ) %>% 
  arrange(desc(admixed_percentage) )


rep_strains_prep <- rep_strains_best_k_long_admix_pops %>%
  group_by(geo, cluster) %>%
  summarise(total_frac = sum(frac_cluster), .groups = 'drop') %>%
  group_by(geo) %>%
  mutate(percent = total_frac / sum(total_frac) * 100) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = levels(pie_best_k_long_admix_pops$cluster)))

stacked_bar_rep_strains <- ggplot(rep_strains_prep, 
                                  aes(x = geo, y = percent, fill = cluster, order = cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = ancestry.colours) +
  labs(
    x = NULL,
    y = "Percentage (%)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none"
  )+ 
  labs(title = "Non-admixed isotypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 10))

stacked_bar_rep_strains



library(ggpubr)

rep_strains_combined_plot <- ggpubr::ggarrange(rep_strains_best_k_admix_plots, 
                                               stacked_bar_rep_strains,
                                   ncol = 1, nrow = 2, 
                                   heights = c(1, 1))
rep_strains_combined_plot

library(dplyr)
library(tidyr)
library(ggplot2)

count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
  count(geo, cluster, name = "n_samples")

p_heatmap_rep_strains<-ggplot(count_df_heatmap_rep_strains, aes(x = cluster, y = geo, fill = n_samples)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    x = "Cluster",
    y = "Geo",
    fill = "Number of Samples",
    title = "Samples per Cluster by Geo"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank()
  )

p_heatmap_rep_strains







library(dplyr)
library(tidyr)
library(ggplot2)

count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
  count(geo, cluster, name = "n_samples") %>%
  mutate(
    geo = factor(geo, levels = sort(unique(geo), decreasing = TRUE))
  )

p_heatmap_rep_strains <- ggplot(count_df_heatmap_rep_strains,
                                aes(x = cluster, y = geo, fill = n_samples)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(n_samples > 0, as.character(n_samples), "")),
            size = 3, color = "black") +
  scale_fill_gradient(low = "lightblue", high = "#1cb6ff") +
  labs(
    x     = NULL,
    y     = NULL,
    fill  = "Number of Samples",
    title = "Non-admixed isotypes per subpoulation by geo"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    plot.title.position   = "plot",      
    legend.position = "bottom",
    plot.title = element_text(
      hjust = 0.65,
      face = "bold",
      size = 10 
    )
    
  )

print(p_heatmap_rep_strains)








library(ggplot2)
library(cowplot)

second_row <- plot_grid(
  stacked_bar, 
  stacked_bar_rep_strains, 
  p_heatmap_rep_strains,
  ncol         = 3,
  rel_widths   = c(1, 1, 2),  
  labels      = c("b","c","d")
)

final_plot <- plot_grid(
  best_k_admix_plots,
  second_row,
  ncol        = 1,
  rel_heights = c(1, 1),
  labels      = c("a","")
)

print(final_plot)
ggsave(paste0("raw_combined_figure_non_cosmopolitan_replicate_",which_replicate,".pdf"), final_plot, width = 10, height = 8, useDingbats = FALSE)



}



