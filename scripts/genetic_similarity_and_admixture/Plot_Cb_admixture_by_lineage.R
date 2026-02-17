
rm(list=ls())



library(pophelper)
library(tidyverse)
library(ggthemes)
library(dplyr)

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





isotype_by_lineage_from_Nic_raw<-read.table("../../processed_data/genetic_similarity_and_admixutre/isotype_byLineage_GeoLocAdmCol_20250909.tsv",
                                            header = TRUE,
                                            sep = '\t',
                                            comment.char  = "")


isotype_by_lineage<-isotype_by_lineage_from_Nic_raw %>% 
  select(isotype,Lineage) %>% 
  filter(isotype %in% sample_names)



best_k_value<-22

which_replicate<-c(1:10)
final_plots <- list()


for (which_replicate in which_replicate) {
  

best_k_qfile<-read.table(paste0("../../processed_data/genetic_similarity_and_admixutre/K22_Processed_Ancestry_replicate_",which_replicate,".tsv"),
                          header = TRUE)

best_k_long_admix_pops <- best_k_qfile %>%
  dplyr::mutate(samples = sample_names) %>%
  tidyr::gather(cluster, frac_cluster, -samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>%
  left_join(isotype_geo_info_ordered, by = c("samples" = "isotype")) %>% 
  left_join(isotype_by_lineage, by = c("samples" = "isotype"))


write.csv(best_k_long_admix_pops,
          paste0("../../processed_data/genetic_similarity_and_admixutre/K22_best_k_long_admix_pops_replicate_",which_replicate,".csv")
)
  

best_k_plot_order <- best_k_long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))


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
        # panel.border=element_rect(linewidth = 0.2),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  coord_cartesian(expand = FALSE) + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))+
  guides(fill = guide_legend(nrow = 2)) + 
  facet_grid(.~ Lineage, scales = "free_x", space = "free_x")+
  theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                    margin = margin(b = 1, t = 3)),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "lines")
        )
  


best_k_admix_plots





pie_best_k_long_admix_pops <- best_k_long_admix_pops %>%
  filter(frac_cluster != "0.000010") %>% 
  group_by(Lineage, cluster) %>%
  summarise(total_frac = sum(frac_cluster), .groups = 'drop') %>%
  group_by(Lineage) %>%
  mutate(percent = total_frac / sum(total_frac) * 100) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = sort(unique(cluster))))

stacked_bar <- ggplot(pie_best_k_long_admix_pops, 
                      aes(x = Lineage, y = percent, fill = cluster, order = cluster)) +
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















rep_strains_best_k_long_admix_pops<-best_k_long_admix_pops %>% 
  filter(max_frac > 0.999)

rep_strains_best_k_plot_order <- rep_strains_best_k_long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))



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
  facet_grid(.~ Lineage, scales = "free_x", space = "free_x")+
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






rep_strains_prep <- rep_strains_best_k_long_admix_pops %>%
  group_by(Lineage, cluster) %>%
  summarise(total_frac = sum(frac_cluster), .groups = 'drop') %>%
  group_by(Lineage) %>%
  mutate(percent = total_frac / sum(total_frac) * 100) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = levels(pie_best_k_long_admix_pops$cluster)))

stacked_bar_rep_strains <- ggplot(rep_strains_prep, 
                                  aes(x = Lineage, y = percent, fill = cluster, order = cluster)) +
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
  count(Lineage, cluster, name = "n_samples")

p_heatmap_rep_strains<-ggplot(count_df_heatmap_rep_strains, aes(x = cluster, y = Lineage, fill = n_samples)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    x = "Cluster",
    y = "Lineage",
    fill = "Number of Samples",
    title = "Samples per Cluster by Lineage"
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
  count(Lineage, cluster, name = "n_samples") %>%
  mutate(
    geo = factor(Lineage, levels = sort(unique(Lineage), decreasing = TRUE))
  )

p_heatmap_rep_strains <- ggplot(count_df_heatmap_rep_strains,
                                aes(x = cluster, y = Lineage, fill = n_samples)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(n_samples > 0, as.character(n_samples), "")),
            size = 3, color = "black") +
  scale_fill_gradient(low = "lightblue", high = "#1cb6ff") +
  labs(
    x     = NULL,
    y     = NULL,
    fill  = "Number of Samples",
    title = "Non-admixed isotypes per subpoulation by lineage"
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
  # align        = "v",
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
if (which_replicate == 5) {
ggsave(paste0("../../figures/SF6_raw_all_715_isotypes_by_lineage_replicate",which_replicate,".pdf"), 
       final_plot, 
       width = 10, 
       height = 8, 
       useDingbats = FALSE)
}













best_k_qfile<-read.table(paste0("../../processed_data/genetic_similarity_and_admixutre/K22_Processed_Ancestry_replicate_",which_replicate,".tsv"),
                         header = TRUE)

best_k_long_admix_pops <- best_k_qfile %>%
  dplyr::mutate(samples = sample_names) %>%
  tidyr::gather(cluster, frac_cluster, -samples) %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(max_frac = max(frac_cluster)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>%
  left_join(isotype_geo_info_ordered, by = c("samples" = "isotype")) %>% 
  left_join(isotype_by_lineage, by = c("samples" = "isotype"))



best_k_long_admix_pops<-best_k_long_admix_pops %>% 
  filter(Lineage %in% c("TD2", "TD3", "ID", "NWD", "Hubei"))



best_k_plot_order <- best_k_long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))

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
                                   size = 3,
                                   margin = margin(t = 0)),    
        axis.text.y=element_blank(),
        axis.title.y = element_text(size = 8,angle = 90, vjust = .5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(linewidth = 0.1),
        axis.ticks.length = unit(0.02, "cm"),
        axis.title.x=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        # panel.border=element_rect(linewidth = 0.2),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  coord_cartesian(expand = FALSE) + 
  theme(legend.position = "right",
        legend.key.size = unit(0.2, "cm"),
  legend.text = element_text(size = 5, margin = margin(l = 0, r = 2)),        # 缩小文字
legend.title = element_text(size = 6), 
legend.key.height = unit(0.25, "cm"),
legend.key.width  = unit(0.25, "cm"), 
legend.spacing.x = unit(0.01, "cm"), 
legend.spacing.y = unit(0, "cm"), 
legend.margin = margin(0,0,0,0)   )+ 
  
  guides(fill = guide_legend(ncol = 3)) +
  facet_grid(.~ Lineage, scales = "free_x", space = "free_x")+
  theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                    margin = margin(b = 1, t = 3)),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")
  )



best_k_admix_plots



rep_strains_best_k_long_admix_pops<-best_k_long_admix_pops %>%
  filter(max_frac > 0.999)

rep_strains_best_k_plot_order <- rep_strains_best_k_long_admix_pops %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::arrange(cluster, -max_frac) %>%
  dplyr::mutate(samples = factor(samples, levels = unique(samples)))



non_admixed_isotype_list<-rep_strains_best_k_plot_order %>% 
  select(-frac_cluster,-max_frac)





library(dplyr)
library(tidyr)
library(ggplot2)

count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
  count(Lineage, cluster, name = "n_samples")


library(dplyr)
library(tidyr)
library(ggplot2)

count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
  count(Lineage, cluster, name = "n_samples") %>%
  mutate(
    geo = factor(Lineage, levels = sort(unique(Lineage), decreasing = TRUE))
  )

p_heatmap_rep_strains <- ggplot(count_df_heatmap_rep_strains,
                                aes(x = cluster, y = Lineage, fill = n_samples)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(n_samples > 0, as.character(n_samples), "")),
            size = 3, color = "black") +
  scale_fill_gradient(low = "lightblue", high = "#1cb6ff") +
  labs(
    x     = NULL,
    y     = NULL,
    fill  = "Number of Samples",
    title = "Non-admixed isotypes per subpoulation by lineage"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    plot.title.position   = "plot",      
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.65,
      face = "bold",
      size = 7 
    ),
    legend.direction = "vertical",
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5),  
    legend.key.width  = unit(0.3, "cm"), 
    legend.key.height = unit(0.35, "cm"), 
    legend.spacing.x  = unit(0.05, "cm") 
  )

print(p_heatmap_rep_strains)




library(ggplot2)
library(cowplot)

final_plot <- plot_grid(
  best_k_admix_plots,
  NULL, 
  p_heatmap_rep_strains,
  ncol = 3,
  rel_widths = c(1, 0.05, 1) 
)

print(final_plot)
# ggsave(paste0("raw_all_715_isotypes_by_lineage_replicate_",which_replicate,"_only_5_lineage.pdf"), final_plot, width = 10, height = 1.4, useDingbats = FALSE)


final_plots[[ paste0("replicate_", which_replicate) ]] <- final_plot

}



final_plots

combined_only_5_lineage <- plot_grid(plotlist = final_plots, ncol = 2, align = "v")
combined_only_5_lineage

ggsave("../../figures/SF7_admixture_only_5_lineage.pdf", 
       combined_only_5_lineage, 
       width = 10, height = 7)

