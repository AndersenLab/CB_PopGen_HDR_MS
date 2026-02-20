rm(list=ls())

library(readxl)
library(geosphere)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)
library(reshape2)
library(ggbreak)
library(patchwork)

Cb_raw_csv<-read.csv("../../supplementary_data/SD1_Cb_2018_strains_data.csv")

#### unique isotypes
iso_unique <- Cb_raw_csv %>%
  dplyr::select(latitude,strain,isotype) %>% 
  dplyr::group_by(isotype) %>%
  dplyr::summarise(n_strains = n_distinct(strain)) %>%
  dplyr::filter(n_strains == 2) %>% 
  dplyr::left_join(Cb_raw_csv, by = c("isotype"))

nrow(iso_unique)


plot_isotype_dist <- function(meta_sheet, condordance_groups, title_text,
                              bins = 20) {
  
  df <- meta_sheet %>%
    dplyr::filter(strain %in% condordance_groups$strain) %>%
    dplyr::select(-isotype) %>%
    dplyr::left_join(condordance_groups, by = "strain")
  
  coords <- df %>%
    dplyr::select(strain, longitude, latitude) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
    dplyr::mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- geosphere::distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- reshape2::melt(dist_mat, varnames = c("strain1", "strain2"), value.name = "geo_m") %>%
    dplyr::filter(!is.na(geo_m)) %>%
    dplyr::mutate(geo_distance = geo_m / 1000) %>%
    dplyr::select(strain1, strain2, geo_distance)
  
  isotype_pairs <- df %>%
    dplyr::select(strain, isotype) %>%
    dplyr::filter(!is.na(isotype)) %>%
    dplyr::group_by(isotype) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::summarise(pairs = list(apply(combn(strain, 2), 2, paste, collapse = "_"))) %>%
    tidyr::unnest(pairs) %>%
    dplyr::pull(pairs)
  
  all_dist <- all_dist %>%
    dplyr::mutate(pair1 = paste(strain1, strain2, sep = "_"),
           pair2 = paste(strain2, strain1, sep = "_"))
  dist_iso <- all_dist %>%
    dplyr::filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>%
    dplyr::pull(geo_distance)
  
  p <- ggplot(data.frame(dist = dist_iso), aes(x = dist)) +
    geom_histogram(bins = bins, fill = "#53886C", color = "black") +
    stat_bin(bins = bins, geom = "text",
             aes(label = after_stat(count)), vjust = -0.5, size = 3) +
    labs(title = title_text,
         x = "Distance (km)",
         y = "Frequency") +
    theme_classic()
  
  
  dist_out <- all_dist %>%
    dplyr::filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>% 
    dplyr::select(strain1, strain2, geo_distance) %>% 
    dplyr::arrange(dplyr::desc(geo_distance))
  
  return(list(plot=p,
              strain_list=dist_out))
}


# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE) %>% 
  dplyr::filter(!(isotype %in% c("MY681","ECA1146","JU356","ECA1503")))

Cb_list <- plot_isotype_dist(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = NULL
)

### Add isotype info
Cb_condordance_groups_tmp<-Cb_condordance_groups_raw %>% 
  dplyr::select(strain,isotype) %>% 
  dplyr::rename(strain1=strain)

Cb_geo_distance_within_same_isotype<- Cb_list$strain_list %>% 
  dplyr::left_join(Cb_condordance_groups_tmp, by =c("strain1") )
  
write.csv(Cb_geo_distance_within_same_isotype,
          "../../processed_data/Geo_info/Cb_geo_distance_within_same_isotype.csv",
          quote = FALSE,
          row.names = FALSE) 


Cb_geo_distance_within_same_isotype %>% 
  dplyr::arrange(desc(geo_distance)) %>% 
  dplyr::select(isotype) %>% 
  dplyr::pull() %>% 
  unique() %>% 
  head()


######### How many isotypes have strains greater than 100 km from each other? 
Cb_geo_distance_within_same_isotype<-read.csv("../../processed_data/Geo_info/Cb_geo_distance_within_same_isotype.csv")

Cb_df_geo_over100km<-Cb_geo_distance_within_same_isotype %>% 
  dplyr::filter(geo_distance >100)

Cb_100km_isotype<-unique(Cb_df_geo_over100km$isotype)
n_Cb_100km_isotype<-length(Cb_100km_isotype)


######## add geo of isotypes > 100 km
indep_info_geo_for_each_strain<-read.csv("../../processed_data/Geo_info/indep_info_geo_for_each_strain.csv")

Cb_over_100km_df<-indep_info_geo_for_each_strain %>% 
  dplyr::filter(isotype %in% Cb_100km_isotype) %>% 
  na.omit()
length(Cb_over_100km_df$isotype %>% unique())
## 38


#### Setting 2000 km as the threshold for defining "Cosmopolitan" samples
Cb_df_geo_over2000km<-Cb_geo_distance_within_same_isotype %>% 
  dplyr::filter(geo_distance >2000)

length(unique(Cb_df_geo_over2000km$isotype))
#17

Cb_Cosmopolitan_isotype<-unique(Cb_df_geo_over2000km$isotype)

write.table(Cb_Cosmopolitan_isotype,
            "../../processed_data/Geo_info/Cb_Cosmopolitan_isotype.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE) 

plot_fs1 <- function(meta_sheet, condordance_groups, title_text,
                     bins = 22) {
  
  df <- meta_sheet %>%
    dplyr::filter(strain %in% condordance_groups$strain) %>%
    dplyr::select(-isotype) %>%
    dplyr::left_join(condordance_groups, by = "strain")
  
  coords <- df %>%
    dplyr::select(strain, longitude, latitude) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
    dplyr::mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- geosphere::distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- reshape2::melt(dist_mat,
                             varnames = c("strain1", "strain2"),
                             value.name = "geo_m") %>%
    dplyr::filter(!is.na(geo_m)) %>%
    dplyr::mutate(geo_distance = geo_m / 1000) %>%
    dplyr::select(strain1, strain2, geo_distance)
  
  isotype_pairs <- df %>%
    dplyr::select(strain, isotype) %>%
    dplyr::filter(!is.na(isotype)) %>%
    dplyr::group_by(isotype) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::summarise(pairs = list(apply(combn(strain, 2), 2, paste, collapse = "_"))) %>%
    tidyr::unnest(pairs) %>%
    dplyr::pull(pairs)
  
  all_dist <- all_dist %>%
    dplyr::mutate(pair1 = paste(strain1, strain2, sep = "_"),
           pair2 = paste(strain2, strain1, sep = "_"))
  
  dist_iso <- all_dist %>%
    dplyr::filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>%
    dplyr::pull(geo_distance)
  
  
  mid_breaks <- seq(100, 5000, length.out = bins - 1)
  breaks <- c(-Inf, mid_breaks, Inf)
  labels <- c(
    "<100",
    paste0(round(mid_breaks[-length(mid_breaks)]), "-", round(mid_breaks[-1])),
    ">5000"
  )
  binned <- cut(dist_iso,
                breaks = breaks,
                labels = labels,
                include.lowest = TRUE,
                right = FALSE)
  
  df_bin <- data.frame(bin = binned) %>%
    count(bin) %>%
    complete(bin = factor(labels, levels = labels), fill = list(n = 0))
  
  p <- ggplot(df_bin, aes(x = bin, y = n)) +
    geom_col(fill = "grey30", color = NA) +
    geom_text(aes(label = n), vjust = -0.5, size = 1.8) +
    labs(title = title_text,
         x = "Distance (km)",
         y = "Frequency") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
          axis.text.y = element_text(size = 6))+
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) 
  
  
  dist_out <- all_dist %>%
    dplyr::filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>%
    dplyr::select(strain1, strain2, geo_distance) %>%
    dplyr::arrange(dplyr::desc(geo_distance))
  

  
  
  p2 <- p +
    scale_y_break(c(6200, 14500),
                  expand = c(0, 100))+
    scale_y_continuous(limits = c(0, 16700),
                       expand = c(0, 0))+
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.line.y.right  = element_blank()
    )+
    theme(
      plot.title      = element_text(size = 8, hjust = 0.5),
      axis.title.x    = element_text(size = 6,
                                     margin = margin(t = -6)),
      axis.title.y    = element_text(size = 6,
                                     margin = margin(r = -6)),
    )

  return(list(plot = p,
              plot_break = p2,
              strain_list = dist_out))
  
}


# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE) %>% 
  filter(!(isotype %in% c("MY681","ECA1146","JU356","ECA1503")))

Cb_fs1 <- plot_fs1(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = NULL
)


#### Max Haversine distance between strains within the same isotype
plot_max_fs1 <- function(meta_sheet, condordance_groups, title_text, bins = 22) {
  df <- meta_sheet %>%
    dplyr::filter(strain %in% condordance_groups$strain) %>%
    dplyr::select(-isotype) %>%
    dplyr::left_join(condordance_groups, by = "strain") %>%
    dplyr::filter(!is.na(isotype))
  
  coords <- df %>%
    dplyr::select(strain, longitude, latitude) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
    dplyr::mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- melt(dist_mat,
                   varnames = c("strain1", "strain2"),
                   value.name = "geo_m") %>%
    dplyr::filter(!is.na(geo_m)) %>%
    dplyr::mutate(geo_distance = geo_m / 1000) %>%
    dplyr::select(strain1, strain2, geo_distance)
 
  iso_info <- df %>% select(strain, isotype)
  dist_iso <- all_dist %>%
    dplyr::left_join(iso_info, by = c("strain1" = "strain")) %>%
    dplyr::rename(isotype1 = isotype) %>%
    dplyr::left_join(iso_info, by = c("strain2" = "strain")) %>%
    dplyr::rename(isotype2 = isotype) %>%
    dplyr::filter(!is.na(isotype1), !is.na(isotype2), isotype1 == isotype2) %>%
    dplyr::group_by(isotype = isotype1) %>%
    dplyr::summarise(max_distance = max(geo_distance)) %>%
    dplyr::ungroup()
  
  dist_max <- dist_iso$max_distance
  
  mid_breaks <- seq(100, 5000, length.out = bins - 1)
  breaks <- c(-Inf, mid_breaks, Inf)
  labels <- c(
    "<100",
    paste0(round(mid_breaks[-length(mid_breaks)]), "-", round(mid_breaks[-1])),
    ">5000"
  )
  binned <- cut(dist_max,
                breaks = breaks,
                labels = labels,
                include.lowest = TRUE,
                right = FALSE)
  
  df_bin <- tibble(bin = binned) %>%
    count(bin) %>%
    complete(bin = factor(labels, levels = labels), fill = list(n = 0))
  
  p <- ggplot(df_bin, aes(x = bin, y = n)) +
    geom_col(fill = "#53886C", color = NA) +
    geom_text(aes(label = n), vjust = -0.5, size = 1.8) +
    labs(title = title_text,
         x = "Max distance (km)",
         y = "Number of isotypes") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))+
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) 
  
  
  p2 <- p +
    scale_y_break(c(18,215),
                  expand = c(0, 1)) +
    scale_y_continuous(limits = c(0, 220), expand = c(0, 0)) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.line.y.right  = element_blank()
    )+
    theme(
      plot.title      = element_text(size = 8, hjust = 0.5),
      axis.title.x    = element_text(size = 6,
                                     margin = margin(t = -6)),
      axis.title.y    = element_text(size = 6,
                                     margin = margin(r = -6)),
    )
  
  return(list(
    plot     = p,
    plot_break = p2,
    max_distances = dist_iso
  ))
}


# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE) %>% 
  filter(!(isotype %in% c("MY681","ECA1146","JU356","ECA1503")))

Cb_fs1_max <- plot_max_fs1(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = NULL
)


#### With annotation Max Haversine distance between strains within the same isotype
plot_max_fs1_annot <- function(meta_sheet, condordance_groups, title_text) {
  df <- meta_sheet %>%
    dplyr::filter(strain %in% condordance_groups$strain) %>%
    dplyr::select(-isotype) %>%
    dplyr::left_join(condordance_groups, by = "strain") %>%
    dplyr::filter(!is.na(isotype))
  
  coords <- df %>%
    dplyr::select(strain, longitude, latitude) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
    dplyr::mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- melt(dist_mat,
                   varnames = c("strain1", "strain2"),
                   value.name = "geo_m") %>%
    dplyr::filter(!is.na(geo_m)) %>%
    dplyr::mutate(geo_distance = geo_m / 1000) %>%
    dplyr::select(strain1, strain2, geo_distance)
  
  iso_info <- df %>% select(strain, isotype)
  dist_iso <- all_dist %>%
    dplyr::left_join(iso_info, by = c("strain1" = "strain")) %>%
    dplyr::rename(isotype1 = isotype) %>%
    dplyr::left_join(iso_info, by = c("strain2" = "strain")) %>%
    dplyr::rename(isotype2 = isotype) %>%
    dplyr::filter(isotype1 == isotype2) %>%
    dplyr::group_by(isotype = isotype1) %>%
    dplyr::summarise(max_distance = max(geo_distance)) %>%
    dplyr::ungroup()
  
  breaks <- c(-Inf, 1, 10, 50, 100, 500, 2000, Inf)
  labels <- c("< 1", "1-10", "10-50", "50-100", "100-500", "500-2000", "> 2000")
  
  df_bin <- dist_iso %>%
    dplyr::mutate(
      bin = cut(max_distance,
                breaks = breaks,
                labels = labels,
                include.lowest = TRUE,
                right = FALSE),
      category = ifelse(max_distance > 2000, "Cosmopolitan", "Non-Cosmopolitan")
    ) %>%
    count(bin, category) %>%
    complete(
      bin = factor(labels, levels = labels),
      category = c("Cosmopolitan", "Non-Cosmopolitan"),
      fill = list(n = 0)
    )
  
  total_bin <- df_bin %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(total = sum(n))
  
  p <- ggplot(df_bin, aes(x = bin, y = n, fill = category)) +
    geom_col() +
    geom_text(
      data = total_bin,
      aes(x = bin, y = total, label = total),
      vjust = -0.3, size = 2.5, inherit.aes = FALSE
    ) +
    scale_fill_manual(
      breaks = c("Non-Cosmopolitan", "Cosmopolitan"),
      values = c("Cosmopolitan" = "#53886C", "Non-Cosmopolitan" = "grey"),
      guide = guide_legend(title = NULL)
    ) +
    labs(
      title = title_text,
      x = "Max distance (km)",
      y = "Number of isotypes"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      legend.position = c(0.8, 0.8),
      legend.key.size = unit(0.8, "lines"),
      legend.text    = element_text(size = 5),
      plot.title      = element_text(size = 8, hjust = 0.5),
      axis.title.x    = element_text(size = 6),
      axis.title.y    = element_text(size = 6)
    ) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8)))+
    scale_y_continuous(limits = c(0, 190), expand = c(0, 0))
    
  
  p2 <- p +
    scale_y_break(c(35, 170), expand = c(0, 1)) +
    scale_y_continuous(limits = c(0, 190), expand = c(0, 0)) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      legend.position = "none",
      plot.title      = element_text(size = 8, hjust = 0.5),
      axis.title.x    = element_text(size = 6, margin = margin(t = -6)),
      axis.title.y    = element_text(size = 6, margin = margin(r = -6))
    )
  
  return(list(
    plot         = p,
    plot_break   = p2,
    max_distances = dist_iso
  ))
}


# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE) %>% 
  dplyr::filter(!(isotype %in% c("MY681","ECA1146","JU356","ECA1503")))


Cb_fs1_max_annot <- plot_max_fs1_annot(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = NULL
)


plot_strain_counts_per_bin <- function(meta_sheet, concordance_groups, title_text) {
  df <- meta_sheet %>%
    dplyr::inner_join(
      concordance_groups %>% select(strain, isotype),
      by = c("strain","isotype")
    )
  
  coords <- df %>%
    dplyr::select(strain, longitude, latitude) %>%
    dplyr::filter(!is.na(longitude), !is.na(latitude)) %>%
    dplyr::mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    ) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- melt(dist_mat,
                   varnames = c("strain1", "strain2"),
                   value.name = "geo_m") %>%
    dplyr::filter(!is.na(geo_m)) %>%
    dplyr::mutate(geo_distance = geo_m / 1000) %>%
    dplyr::select(strain1, strain2, geo_distance)
  
  iso_info <- df %>% dplyr::select(strain, isotype)
  dist_iso <- all_dist %>%
    dplyr::left_join(iso_info, by = c("strain1" = "strain")) %>%
    dplyr::rename(isotype1 = isotype) %>%
    dplyr::left_join(iso_info, by = c("strain2" = "strain")) %>%
    dplyr::rename(isotype2 = isotype) %>%
    dplyr::filter(isotype1 == isotype2) %>%
    dplyr::group_by(isotype = isotype1) %>%
    dplyr::summarise(max_distance = max(geo_distance)) %>%
    dplyr::ungroup()
  
  breaks <- c(-Inf, 1, 10, 50, 100, 500, 2000, Inf)
  labels <- c("< 1", "1-10", "10-50", "50-100", "100-500", "500-2000", "> 2000")
  
  strain_bins <- df %>%
    dplyr::select(strain, isotype) %>%
    dplyr::left_join(dist_iso, by = "isotype") %>%
    dplyr::mutate(bin = cut(max_distance,
                     breaks = breaks,
                     labels = labels,
                     include.lowest = TRUE,
                     right = FALSE))  %>% 
    na.omit() %>%
    dplyr::mutate(category = ifelse(max_distance > 100, "Cosmopolitan", "Non-Cosmopolitan")
           )
  
  strain_count <- strain_bins %>%
    count(bin, name = "n_strains") %>%
    complete(bin = factor(labels, levels = labels),
             fill = list(n_strains = 0)) %>% 
    dplyr::mutate(category = ifelse(bin %in% c(
      "> 2000"),
                             "Cosmopolitan", "Non-Cosmopolitan")) %>% 
    dplyr::mutate(category = factor(category,
                             levels = c("Non-Cosmopolitan","Cosmopolitan")))
  

  
    p_strains <- ggplot(strain_count, aes(x = bin, 
                                          y = n_strains,
                                          fill = category)) +
    geom_col() +
    geom_text(aes(label = n_strains),
              vjust = -0.3, size = 2.5) +
      scale_fill_manual(
        values = c(
          "Non-Cosmopolitan" = "grey80",
          "Cosmopolitan"     = "#53886C"
        ),guide = guide_legend(title = NULL)
      )+
    labs(
      title = title_text,
      x     = "Max distance (km)",
      y     = "Number of strains"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title  = element_text(size = 8, hjust = 0.5),
      axis.title  = element_text(size = 6),
      legend.position = c(0.8, 0.8),
      legend.key.size = unit(0.8, "lines"),
      legend.text    = element_text(size = 5)
    ) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) +
    scale_y_continuous(limits = c(0, 730),
                       expand = c(0, 0))
  
  return(p_strains)
}


Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE) %>% 
  dplyr::filter(!(isotype %in% c("MY681","ECA1146","JU356","ECA1503")))


Cb_fs1_max_annot_strains <- plot_strain_counts_per_bin(
  meta_sheet = Cb_raw_csv,
  concordance_groups = Cb_condordance_groups_raw,
  title_text = NULL
)


p1 <- Cb_fs1_max_annot$plot
p2 <- Cb_fs1_max_annot_strains

combined <- ggpubr::ggarrange(ggpubr::ggarrange(p1,
                                                p2,
                                         ncol = 2, 
                                         labels = c("c","d"),
                                         widths = c(0.5, 0.5)))

saveRDS(combined,"../../processed_data/genetic_similarity_and_admixutre/geo_distance_plots.RData")


##### plot map of these 18 isotypes ###
world <- map_data('world')
world_map <- world[world$region != "Antarctica",] # intercourse antarctica

Cb_Cosmopolitan_isotype_info<-Cb_raw_csv %>% 
  dplyr::select(strain,isotype,latitude,longitude) %>% 
  dplyr::filter(isotype %in% Cb_Cosmopolitan_isotype) 
length(unique(Cb_Cosmopolitan_isotype_info$isotype))
#17


p_map_Cosmopolitan <- ggplot(Cb_Cosmopolitan_isotype_info, aes(x = longitude, y = latitude)) +
  geom_map(data = world_map, map = world_map,
           aes(map_id = region, x = long, y = lat),
           fill = "#F7F7F7", color = "black", linewidth = 0.03) +
  geom_point(shape = 16, size = 1, color = "#53886C", alpha = 0.8) +
  lims(x = c(-180, 191), y = c(-58, 84)) +
  facet_wrap(~ isotype, ncol = 6) +
  theme_void() +
  theme(strip.text = element_text(size = 6))


p_map_Cosmopolitan_3_col <- ggplot(Cb_Cosmopolitan_isotype_info, aes(x = longitude, y = latitude)) +
  geom_map(data = world_map, map = world_map,
           aes(map_id = region, x = long, y = lat),
           color="gray95", fill="gray80", linewidth=0.05)+
  geom_point(shape = 16, size = 2, color = "red", alpha = 0.5) +
  lims(x = c(-180, 191), y = c(-58, 84)) +
  facet_wrap(~ isotype, ncol = 3) +
  theme_void() +
  theme(strip.text = element_text(size = 8))


label_positions <- Cb_Cosmopolitan_isotype_info %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(longitude = -170,
            latitude = -50)

p_map_Cosmopolitan_4_col <- ggplot(Cb_Cosmopolitan_isotype_info, aes(x = longitude, y = latitude)) +
  geom_map(data = world_map, map = world_map,
           aes(map_id = region, x = long, y = lat),
           color="gray95", fill="gray80", linewidth=0.02)+
  geom_point(shape = 16, size = 2, color = "red", alpha = 0.5) +
  lims(x = c(-180, 191), y = c(-58, 84)) +
  facet_wrap(~ isotype, ncol = 4) +
  geom_text(data = label_positions, 
            aes(x = longitude, y = latitude, label = isotype),
            hjust = 0, vjust = 0, size = 3, fontface = "bold") +
  theme_void() +
  theme(strip.text = element_blank())+
  theme(panel.spacing = unit(0.05, "inches"))


isotype_counts <- Cb_Cosmopolitan_isotype_info %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(n_strains = n())

Cb_Cosmopolitan_strain_frequency_plot<-ggplot(isotype_counts, aes(x = reorder(isotype, -n_strains), y = n_strains)) +
  geom_col(fill = "red",width = 0.7) +
  geom_text(aes(label = n_strains), vjust = -0.5, size = 2.5) +
  labs(
       y = "Number of strains") +
  theme_classic()+
  scale_y_continuous(limits = c(0, 330))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7), 
    plot.title = element_text(hjust = 0.5, size = 6) 
  )


#### merge maps and frequency
library(cowplot)

combined_plot <- ggdraw() +
  draw_plot(p_map_Cosmopolitan_4_col, 0, 0, 1, 1) + 
  draw_plot(Cb_Cosmopolitan_strain_frequency_plot, 0.25, 0, 0.75, 0.2) 

save(combined_plot, file = "../../processed_data/genetic_similarity_and_admixutre/Cosmo_18_maps.RData")
