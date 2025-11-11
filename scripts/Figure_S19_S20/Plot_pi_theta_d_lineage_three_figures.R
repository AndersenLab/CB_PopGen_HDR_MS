rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(cowplot)

source("../utilities.R")

lineage<-readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byLineage_GeoLocAdmCol_20250909.tsv")
table(lineage$Lineage)
# AD        FM     Hubei        ID        KD       NWD 
# 32         4         2         3        10         4 
# Quebec       TD1       TD2       TD3 Temperate        TH 
# 1        27         7         3        29        95 
# Tropical 
# 502 


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



# lineage_display_names[lineage_display_names == "Cb_non_globalized"] <- "Non-globalized"
# lineage_display_names[lineage_display_names == "Asia_div"] <- "Asia div"
# lineage_display_names[lineage_display_names == "Asia_non_div"] <- "Asia non-div"
# lineage_display_names[lineage_display_names == "Australia_div"] <- "Australia div"
# lineage_display_names[lineage_display_names == "Australia_non_div"] <- "Australia non-div"
# lineage_display_names[lineage_display_names == "Central_America"] <- "Central America"
# lineage_display_names[lineage_display_names == "South_America"] <- "South America"
names(lineage_display_names) <- names(lineage_info)




##### read genome domain
genome_domain_raw<-read.table("../../data/05.07.21_cb_subregion.bed",
                            header = TRUE,sep = '\t')

genome_domain <- genome_domain_raw %>%
  # merge sub_lineages
  mutate(
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
  filter(!is.na(category)) %>%
  mutate(category = factor(category, levels = c("Tip", "Arm", "Center"))) %>% 
  rename(chrom = CHROM)













###### pi ##########


# Load and combine all data
all_pi_df <- purrr::map2_dfr(names(lineage_info), lineage_info, function(lineage, base_dir) {
  file_path <- file.path(base_dir, 
                         # lineage, 
                         "chromosome_windows_diversity.csv")
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(stat_type == "pi") %>%
    mutate(lineage = lineage_display_names[[lineage]])
})

# # Normalize position within each chrom and lineage
# all_pi_df <- all_pi_df %>%
#   group_by(lineage, chrom) %>%
#   mutate(x_norm = (x - min(x)) / (max(x) - min(x))) %>%
#   ungroup()








# Order lineages so that the facet grid knows the bottom row
lineage_levels <- unique(lineage_display_names[names(lineage_info)])
all_pi_df$lineage <- factor(all_pi_df$lineage, levels = lineage_levels)

# Create plot
pi_plot <- ggplot(
  # all_pi_df, aes(x = x/1e6, y = stat)
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
      # face = "bold", 
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

ggsave("../../figures/FigureS19_all_lineage_pi_same_scale.pdf", pi_plot, width = 7, height = 7)




# pi_plot_free_y <- ggplot(all_pi_df, aes(x = x_norm, y = stat)) +
#   geom_point(color = "gray30", size = 0.2, alpha = 0.8, shape = 16) +
#   geom_smooth(method = "loess", se = FALSE, span = 0.3, color = "lightgray") +
#   facet_grid(lineage ~ chrom, scales = "free_y") +
#   theme_bw(base_size = 9) +
#   theme(
#     panel.grid = element_blank(),
#     strip.background = element_blank(),
#     strip.text = element_text(
#       # face = "bold", 
#       size = 6),
#     axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 7),
#     axis.ticks.x = element_line(),
#     axis.title = element_text(size = 9),
#     plot.margin = margin(2, 2, 2, 2)
#   ) +
#   labs(
#     x = "Normalized genome position",
#     y = expression("Nucleotide diversity (" * pi * ")")
#   )
# 
# ggsave("all_lineage_pi_free_y.pdf", pi_plot_free_y, width = 7.5, height = 10)

















###### theta ##########

# Load and combine all data
all_theta_df <- purrr::map2_dfr(names(lineage_info), lineage_info, function(lineage, base_dir) {
  file_path <- file.path(base_dir, 
                         # lineage, 
                         "chromosome_windows_diversity.csv")
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(stat_type == "theta") %>%
    mutate(lineage = lineage_display_names[[lineage]])
})

# # Normalize position within each chrom and lineage
# all_theta_df <- all_theta_df %>%
#   group_by(lineage, chrom) %>%
#   mutate(x_norm = (x - min(x)) / (max(x) - min(x))) %>%
#   ungroup()

# Order lineages so that the facet grid knows the bottom row
lineage_levels <- unique(lineage_display_names[names(lineage_info)])
all_theta_df$lineage <- factor(all_theta_df$lineage, levels = lineage_levels)

# Create plot
theta_plot <- ggplot(
  # all_pi_df, aes(x = x/1e6, y = stat)
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
      # face = "bold", 
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

ggsave("../../figures/FigureS20_all_lineage_theta_same_scale.pdf", theta_plot, width = 7, height = 7)






