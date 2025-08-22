rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(cowplot)

source("../utilities.R")



# Define all regions and corresponding base directories
region_info <- list(
  `./` = "../../processed_data/pi_theta_d_1kb",
  # cosmopolitan = "../../processed_data/pi_theta_d_geo",
  Cb_non_cosmopolitan = "../../processed_data/pi_theta_d_geo",
  Asia = "../../processed_data/pi_theta_d_geo",
  # Asia_div = "../../processed_data/pi_theta_d_separate_Au_As",
  # Asia_non_div = "../../processed_data/pi_theta_d_separate_Au_As",
  Australia = "../../processed_data/pi_theta_d_geo",
  # Australia_div = "../../processed_data/pi_theta_d_separate_Au_As",
  # Australia_non_div = "../../processed_data/pi_theta_d_separate_Au_As",
  Caribbean = "../../processed_data/pi_theta_d_geo",
  Central_America = "../../processed_data/pi_theta_d_geo",
  Hawaii = "../../processed_data/pi_theta_d_geo",
  # North_America = "../../processed_data/pi_theta_d_geo",
  Pacific = "../../processed_data/pi_theta_d_geo",
  South_America = "../../processed_data/pi_theta_d_geo",
  Taiwan = "../../processed_data/pi_theta_d_geo"
)



### sample size
geo_raw<-read.csv("../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv")

sample_size<-geo_raw %>% 
  group_by(geo) %>% 
  summarise(n=n()) %>% 
  ungroup()

sample_size
# 1 Africa             23
# 2 Asia               60
# 3 Australia          57
# 4 Caribbean          41
# 5 Central America   121
# 6 Europe             12
# 7 Hawaii             80
# 8 North America      18
# 9 Pacific            50
# 10 South America      62
# 11 Taiwan            153
# 12 cosmopolitan       38
# 13 unknown             4


## All
nrow(geo_raw)
# 719

## Non-cosmopolitan
nrow(geo_raw)-sample_size[which(sample_size$geo=="cosmopolitan"),2]
# 681



# Display name mapping
region_display_names <- names(region_info)
region_display_names[region_display_names == "./"] <- "All (719)"
region_display_names[region_display_names == "Cb_non_cosmopolitan"] <- "Non-cosmopolitan (681)"
# region_display_names[region_display_names == "Africa"] <- "Africa (23)"
region_display_names[region_display_names == "Asia"] <- "Asia (60)"
region_display_names[region_display_names == "Australia"] <- "Australia (57)"
region_display_names[region_display_names == "Caribbean"] <- "Caribbean (57)"
region_display_names[region_display_names == "Central_America"] <- "Central America (121)"
region_display_names[region_display_names == "Hawaii"] <- "Hawaii (80)"
region_display_names[region_display_names == "Pacific"] <- "Hawaii (50)"
region_display_names[region_display_names == "South_America"] <- "South America (62)"
region_display_names[region_display_names == "Taiwan"] <- "Taiwan (153)"

names(region_display_names) <- names(region_info)




##### read genome domain
genome_domain_raw<-read.csv("../../data/Cb_bounds_df.csv",
                            header = TRUE)

genome_domain <- genome_domain_raw %>%
  # merge sub_regions
  mutate(
    category = case_when(
      grepl("tip$", sub_region)    ~ "Tip",
      grepl("arm$", sub_region)    ~ "Arm",
      sub_region == "center"       ~ "Center",
      TRUE                         ~ NA_character_
    ),
    # as Mb 
    xmin = start / 1e6,
    xmax = stop  / 1e6
  ) %>%
  filter(!is.na(category)) %>%
  mutate(category = factor(category, levels = c("Tip", "Arm", "Center")))












###### pi ##########


# Load and combine all data
all_pi_df <- purrr::map2_dfr(names(region_info), region_info, function(region, base_dir) {
  file_path <- file.path(base_dir, region, "chromosome_windows_diversity.csv")
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(stat_type == "pi") %>%
    mutate(region = region_display_names[[region]])
})

# # Normalize position within each chrom and region
# all_pi_df <- all_pi_df %>%
#   group_by(region, chrom) %>%
#   mutate(x_norm = (x - min(x)) / (max(x) - min(x))) %>%
#   ungroup()








# Order regions so that the facet grid knows the bottom row
region_levels <- unique(region_display_names[names(region_info)])
all_pi_df$region <- factor(all_pi_df$region, levels = region_levels)

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
              method = "loess", se = FALSE, span = 0.3, color = "lightgray") +
  facet_grid(region ~ chrom, scales = "free_x") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      # face = "bold", 
      size = 5),
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_line(),
    axis.title = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  theme(legend.position = "none")+
  ylim(0, 0.02)+  ##### removed too divergent sites
  labs(
    x = "Physical genome position",
    y = expression("Nucleotide diversity (" * pi * ")")
  )

# Save

# ggsave("Max_0.02_all_region_pi_same_scale.pdf", pi_plot, width = 7, height = 9)
ggsave("Max_0.02_all_region_pi_same_scale.png", pi_plot, width = 7, height = 9, dpi = 300)



















###### theta ##########

# Load and combine all data
all_theta_df <- purrr::map2_dfr(names(region_info), region_info, function(region, base_dir) {
  file_path <- file.path(base_dir, region, "chromosome_windows_diversity.csv")
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(stat_type == "theta") %>%
    mutate(region = region_display_names[[region]])
})

# Order regions so that the facet grid knows the bottom row
region_levels <- unique(region_display_names[names(region_info)])
all_theta_df$region <- factor(all_theta_df$region, levels = region_levels)

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
              method = "loess", se = FALSE, span = 0.3, color = "lightgray") +
  facet_grid(region ~ chrom, scales = "free_x") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      # face = "bold", 
      size = 5),
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_line(),
    axis.title = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  ylim(0, 0.02)+  ##### removed too divergent sites
  theme(legend.position = "none")+
  labs(
    x = "Physical genome position",
    y = expression("Watterson's " * theta * "")
  )

# Save

# ggsave("all_region_theta_same_scale.pdf", theta_plot, width = 7, height = 9)
ggsave("Max_0.02_all_region_theta_same_scale.png", theta_plot, width = 7, height = 9,dpi = 300)






















###### d ##########

# Load and combine all data
all_d_df <- purrr::map2_dfr(names(region_info), region_info, function(region, base_dir) {
  file_path <- file.path(base_dir, region, "chromosome_windows_diversity.csv")
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(stat_type == "d") %>%
    mutate(region = region_display_names[[region]])
})

# Normalize position within each chrom and region
all_d_df <- all_d_df %>%
  group_by(region, chrom) %>%
  mutate(x_norm = (x - min(x)) / (max(x) - min(x))) %>%
  ungroup()

# Order regions so that the facet grid knows the bottom row
region_levels <- unique(region_display_names[names(region_info)])
all_d_df$region <- factor(all_d_df$region, levels = region_levels)

# Create plot
d_plot <- ggplot(
  # all_pi_df, aes(x = x/1e6, y = stat)
) +
  geom_rect(
    data = genome_domain,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  scale_fill_manual(values = genome_domain_colors) +
  geom_point(data=all_d_df, aes(x = x/1e6, y = stat),
             color = "gray30", size = 0.2, alpha = 0.8, shape = 16) +
  geom_smooth(data=all_d_df, aes(x = x/1e6, y = stat),
              method = "loess", se = FALSE, span = 0.3, color = "lightgray") +
  facet_grid(region ~ chrom, scales = "free_x") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      # face = "bold", 
      size = 5),
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_line(),
    axis.title = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  theme(legend.position = "none")+
  labs(
    x = "Physical genome position",
    y = "Tajima's D"
  )

# Save

# ggsave("all_region_d_same_scale.pdf", d_plot, width = 7, height = 9)
ggsave("all_region_d_same_scale.png", d_plot, width = 7, height = 9,dpi = 300)










