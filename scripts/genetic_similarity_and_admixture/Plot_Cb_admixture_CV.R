rm(list=ls())

library(pophelper)
library(tidyverse)
library(ggthemes)
library(dplyr)

source("../utilities.R")

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

desired <- as.character(sample_list$V1)

isotype_geo_info_ordered <- isotype_geo_info %>%
  dplyr::filter(isotype %in% desired) %>%
  dplyr::mutate(isotype = factor(isotype, levels = desired)) %>%
  dplyr::arrange(isotype)

# get list of isotype names 
sample_names <- isotype_geo_info_ordered$isotype
best_k <- data.frame(K.22 = 20:24)
k_values <- best_k$K.22
best_k_value <- as.numeric(sub("K.", "", colnames(best_k)[1]))

which_replicate<-c(1:10)

for (which_replicate in which_replicate) {
  admix_plots <- list()

for(kpops in 1:length(grep(".Q", list.files("../../processed_data/genetic_similarity_and_admixutre/Q_files"), value = T))){
  K <- as.numeric(sub(".Q", "", strsplit(grep(".Q", list.files("../../processed_data/genetic_similarity_and_admixutre/Q_files"), value = T)[kpops], split = "\\_")[[1]][3]))

  if (!(K %in% k_values)) {
    next
  }

  # load Q files
  qfile_name <- grep(pattern = glue::glue("_{K}_\\d+\\.Q$"), value = T, x = list.files("../../processed_data/genetic_similarity_and_admixutre/Q_files/"))
  qfile <- pophelper::readQ(files = paste0("../../processed_data/genetic_similarity_and_admixutre/Q_files/",qfile_name))[[which_replicate]]
  # add pop names
  names_pool <- c(LETTERS, paste0("A", LETTERS))
  colnames(qfile) <- names_pool[1:K]

  qfile <- qfile %>%
    dplyr::mutate(samples = sample_names)

  if (K == best_k_value) {
    write.table(qfile, file = paste0("../../processed_data/genetic_similarity_and_admixutre/K", K, "_Processed_Ancestry_replicate_",which_replicate,".tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
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

final_admixture_plots <- cowplot::plot_grid(
  admixture_plots,
  legend,
  ncol = 1,
  rel_heights = c(4, 0.3)
)

if (which_replicate == 5) {
  ggsave(
    plot = final_admixture_plots,
    filename = paste0("../../figures/SF5_admixture_replicate_", which_replicate, ".pdf"),
    height = 7,
    width = 7,
    useDingbats = FALSE
  )
  }
}

