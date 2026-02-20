rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(cowplot)

source("../utilities.R")

lineage<-readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv")

# Define all lineages and corresponding base directories
lineage_info <- list(
  `./` = "../../processed_data/diversity_and_divergence/pi_theta_d",
  AD ="../../processed_data/diversity_and_divergence/pi_theta_d_rg/AD/",
  KD ="../../processed_data/diversity_and_divergence/pi_theta_d_rg/KD/",
  TD1 = "../../processed_data/diversity_and_divergence/pi_theta_d_rg/TD1",
  Temperate="../../processed_data/diversity_and_divergence/pi_theta_d_rg/Temperate",
  TH="../../processed_data/diversity_and_divergence/pi_theta_d_rg/TH",
  Tropical="../../processed_data/diversity_and_divergence/pi_theta_d_rg/Tropical"
)

# Display name mapping
lineage_display_names <- names(lineage_info)
lineage_display_names[lineage_display_names == "./"] <- "All (715)"
lineage_display_names[lineage_display_names == "AD"] <- "AD (32)"
lineage_display_names[lineage_display_names == "KD"] <- "KD (10)"
lineage_display_names[lineage_display_names == "TD1"] <- "TD1 (27)"
lineage_display_names[lineage_display_names == "Temperate"] <- "Temperate (29)"
lineage_display_names[lineage_display_names == "TH"] <- "TH (95)"
lineage_display_names[lineage_display_names == "Tropical"] <- "Tropical (502)"

names(lineage_display_names) <- names(lineage_info)

##### read genome domain
genome_domain_raw<-read.table("../../data/briggsae_genome_files/05.07.21_cb_subregion.bed",
                            header = TRUE,sep = '\t')

genome_domain <- genome_domain_raw %>%
  # merge sub_lineages
  dplyr::mutate(
    category = case_when(
      grepl("tip$", Location)    ~ "Tip",
      grepl("arm$", Location)    ~ "Arm",
      Location == "center"       ~ "Center",
      TRUE                         ~ NA_character_
    ),
    # as Mb 
    xmin = start / 1e6,
    xmax = end  / 1e6
  ) %>%
  dplyr::filter(!is.na(category)) %>%
  dplyr::mutate(category = factor(category, levels = c("Tip", "Arm", "Center"))) %>% 
  rename(chrom = CHROM)


###### pi ##########

# Load and combine all data
all_pi_df <- purrr::map2_dfr(names(lineage_info), lineage_info, function(lineage, base_dir) {
  file_path <- file.path(base_dir, 
                         "chromosome_windows_diversity.csv")
  readr::read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::filter(stat_type == "pi") %>%
    dplyr::mutate(lineage = lineage_display_names[[lineage]])
})

# Order lineages so that the facet grid knows the bottom row
lineage_levels <- unique(lineage_display_names[names(lineage_info)])
all_pi_df$lineage <- factor(all_pi_df$lineage, levels = lineage_levels)

# Create plot
pi_plot <- ggplot(
  ) +
  geom_rect(
    data = genome_domain,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  scale_fill_manual(values = genome_domain_colors) +
  geom_point(data=all_pi_df, aes(x = x/1e6, y = stat),
             color = "gray30", size = 0.2, alpha = 0.8, shape = 16) +
  geom_smooth(data=all_pi_df, aes(x = x/1e6, y = stat),
              method = "loess", se = FALSE, span = 0.3, color = "lightgray",
              linewidth = 0.5, alpha = 0.2) +
  facet_grid(lineage ~ chrom, scales = "free_x") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 6),
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_line(),
    axis.title = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  theme(legend.position = "none")+
  labs(
    x = "Physical genome position",
    y = expression("Nucleotide diversity (" * pi * ")")
  )

# Save
ggsave("../../figures/SF12_all_lineage_pi_same_scale.pdf", pi_plot, width = 7, height = 7)


###### theta ##########

# Load and combine all data
all_theta_df <- purrr::map2_dfr(names(lineage_info), lineage_info, function(lineage, base_dir) {
  file_path <- file.path(base_dir, 
                         # lineage, 
                         "chromosome_windows_diversity.csv")
  readr::read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::filter(stat_type == "theta") %>%
    dplyr::mutate(lineage = lineage_display_names[[lineage]])
})

# Order lineages so that the facet grid knows the bottom row
lineage_levels <- unique(lineage_display_names[names(lineage_info)])
all_theta_df$lineage <- factor(all_theta_df$lineage, levels = lineage_levels)

# Create plot
theta_plot <- ggplot(
) +
  geom_rect(
    data = genome_domain,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  scale_fill_manual(values = genome_domain_colors) +
  geom_point(data=all_theta_df, aes(x = x/1e6, y = stat),
             color = "gray30", size = 0.2, alpha = 0.8, shape = 16) +
  geom_smooth(data=all_theta_df, aes(x = x/1e6, y = stat),
              method = "loess", se = FALSE, span = 0.3, color = "lightgray",
              linewidth = 0.5, alpha = 0.2) +
  facet_grid(lineage ~ chrom, scales = "free_x") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 6),
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_line(),
    axis.title = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  theme(legend.position = "none")+
  labs(
    x = "Physical genome position",
    y = expression("Watterson's " * theta * "")
  )

# Save
ggsave("../../figures/SF13_all_lineage_theta_same_scale.pdf", theta_plot, width = 7, height = 7)

