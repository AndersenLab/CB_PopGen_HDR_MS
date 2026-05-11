rm(list = ls())

library(geosphere)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)

# ============================================================
# Input data
# ============================================================

Cb_raw_csv <- readr::read_tsv("../../data/master_sheet_cb/Cb_1972_strains_data.tsv") %>%
  filter(!(isotype %in% c("MY681", "ECA1146", "JU356", "ECA1503"))) %>%
  mutate(
    longitude = ifelse(
      strain %in% c("PB857", "PB858", "PB859"),
      "-84.0554",
      longitude
    )
  )

Cb_condordance_groups_raw <- read.table(
  "../../data/briggsae_isotypes/isotype_groups.tsv",
  header = TRUE
) %>%
  filter(!(isotype %in% c("MY681", "ECA1146", "JU356", "ECA1503")))


# ============================================================
# Helper: calculate max geographic distance within each isotype
# ============================================================

calculate_max_distance_by_isotype <- function(meta_sheet, concordance_groups) {
  
  df <- meta_sheet %>%
    inner_join(
      concordance_groups %>% select(strain, isotype),
      by = c("strain", "isotype")
    ) %>%
    filter(!is.na(isotype))
  
  coords <- df %>%
    select(strain, longitude, latitude) %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    ) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- geosphere::distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- reshape2::melt(
    dist_mat,
    varnames = c("strain1", "strain2"),
    value.name = "geo_m"
  ) %>%
    filter(!is.na(geo_m)) %>%
    mutate(geo_distance = geo_m / 1000) %>%
    select(strain1, strain2, geo_distance)
  
  iso_info <- df %>%
    select(strain, isotype)
  
  dist_iso <- all_dist %>%
    left_join(iso_info, by = c("strain1" = "strain")) %>%
    rename(isotype1 = isotype) %>%
    left_join(iso_info, by = c("strain2" = "strain")) %>%
    rename(isotype2 = isotype) %>%
    filter(!is.na(isotype1), !is.na(isotype2), isotype1 == isotype2) %>%
    group_by(isotype = isotype1) %>%
    summarise(max_distance = max(geo_distance), .groups = "drop")
  
  list(
    df = df,
    dist_iso = dist_iso
  )
}


# ============================================================
# Panel c: number of isotypes by max distance bin
# ============================================================

plot_max_distance_isotype_count <- function(dist_iso, title_text = NULL) {
  
  breaks <- c(-Inf, 1, 10, 50, 100, 500, 2000, Inf)
  labels <- c("< 1", "1-10", "10-50", "50-100", "100-500", "500-2000", "> 2000")
  
  df_bin <- dist_iso %>%
    mutate(
      bin = cut(
        max_distance,
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    count(bin, name = "n") %>%
    complete(
      bin = factor(labels, levels = labels),
      fill = list(n = 0)
    )
  
  p <- ggplot(df_bin, aes(x = bin, y = n)) +
    geom_col(fill = "grey") +
    geom_text(aes(label = n), vjust = -0.3, size = 2.5) +
    labs(
      title = title_text,
      x = "Max distance (km)",
      y = "Number of isotypes"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 8, hjust = 0.5),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6)
    ) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) +
    scale_y_continuous(limits = c(0, 190), expand = c(0, 0))
  
  p
}


# ============================================================
# Panel d: number of strains by max distance bin
# ============================================================

plot_max_distance_strain_count <- function(df, dist_iso, title_text = NULL) {
  
  breaks <- c(-Inf, 1, 10, 50, 100, 500, 2000, Inf)
  labels <- c("< 1", "1-10", "10-50", "50-100", "100-500", "500-2000", "> 2000")
  
  strain_count <- df %>%
    select(strain, isotype) %>%
    left_join(dist_iso, by = "isotype") %>%
    mutate(
      bin = cut(
        max_distance,
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    filter(!is.na(bin)) %>%
    count(bin, name = "n_strains") %>%
    complete(
      bin = factor(labels, levels = labels),
      fill = list(n_strains = 0)
    )
  
  p <- ggplot(strain_count, aes(x = bin, y = n_strains)) +
    geom_col(fill = "grey") +
    geom_text(aes(label = n_strains), vjust = -0.3, size = 2.5) +
    labs(
      title = title_text,
      x = "Max distance (km)",
      y = "Number of strains"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 8, hjust = 0.5),
      axis.title = element_text(size = 6)
    ) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) +
    scale_y_continuous(limits = c(0, 730), expand = c(0, 0))
  
  p
}


# ============================================================
# Generate combined plot and save RDS
# ============================================================

geo_distance_data <- calculate_max_distance_by_isotype(
  meta_sheet = Cb_raw_csv,
  concordance_groups = Cb_condordance_groups_raw
)

p1 <- plot_max_distance_isotype_count(
  dist_iso = geo_distance_data$dist_iso,
  title_text = NULL
)

p2 <- plot_max_distance_strain_count(
  df = geo_distance_data$df,
  dist_iso = geo_distance_data$dist_iso,
  title_text = NULL
)

combined <- ggpubr::ggarrange(
  p1,
  p2,
  ncol = 2,
  labels = c("c", "d"),
  widths = c(0.5, 0.5)
)

saveRDS(
  combined,
  "../../processed_data/genetic_similarity_and_admixutre/geo_distance_plots.RData"
)



single_strain_isotypes_2 <- Cb_raw_csv %>%
  select(strain,isotype,latitude) %>% 
  count(isotype, name = "n_strain") %>%
  filter(n_strain == 2)

test_2 <- Cb_raw_csv %>%
  filter(isotype %in% single_strain_isotypes_2$isotype) %>% 
  select(strain,isotype,latitude)


single_strain_isotypes_3 <- Cb_raw_csv %>%
  select(strain,isotype,latitude) %>% 
  count(isotype, name = "n_strain") %>%
  filter(n_strain == 3)

test_3 <- Cb_raw_csv %>%
  filter(isotype %in% single_strain_isotypes_3$isotype) %>% 
  select(strain,isotype,latitude)


