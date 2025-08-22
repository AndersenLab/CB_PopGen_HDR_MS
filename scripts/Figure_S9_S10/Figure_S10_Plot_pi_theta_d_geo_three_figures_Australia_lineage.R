rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(stringr)


source("../utilities.R")
# —————————————————————————————————————————
# （1）all isotypes of the geo 
# —————————————————————————————————————————
geo_df <- read.csv(
  "../../processed_data/pi_theta_d_geo/Australia/chromosome_windows_diversity.csv",
  header = TRUE
) %>%
  filter(stat_type %in% c("pi", "theta", "d")) %>%
  mutate(
    region    = "All Australia",
    stat_type = recode(stat_type,
                       pi    = "Nucleotide diversity (\u03C0)",
                       theta = "Watterson's \u03B8",
                       d     = "Tajima's D")
  ) %>%
  group_by(chrom) %>%
  mutate(
    x = x / 1e6,
    x_norm = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  ) %>%
  ungroup()
# —————————————————————————————————————————
# （2）lineage function
# —————————————————————————————————————————
read_window_data <- function(region_name,
                             base_dir = "../../processed_data/pi_theta_d_by_lineage/Australia_lineage") {
  win_file <- file.path(base_dir, region_name, "chromosome_windows_diversity.csv")
  read.csv(win_file, header = TRUE, row.names = 1) %>%
    filter(stat_type %in% c("pi", "theta", "d")) %>%
    mutate(
      region    = region_name,
      stat_type = recode(stat_type,
                         pi    = "Nucleotide diversity (\u03C0)",
                         theta = "Watterson's \u03B8",
                         d     = "Tajima's D")
    ) %>%
    group_by(chrom) %>%
    mutate(
      x = x / 1e6,
      x_norm = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    ) %>%
    ungroup()
}

# —————————————————————————————————————————
# （3）read lineages & merge
# —————————————————————————————————————————
regions <- c("AD", "Tropical")
all_diversity <- lapply(regions, read_window_data) %>% bind_rows()

all_diversity <- bind_rows(geo_df, all_diversity) %>%
  mutate(
    region = factor(region,
                    levels = c("All Australia", regions)),
    stat_type = factor(stat_type,
                       levels = c("Nucleotide diversity (\u03C0)",
                                  "Watterson's \u03B8",
                                  "Tajima's D"))
  )

# —————————————————————————————————————————
# （4）read genome domain 
# —————————————————————————————————————————
genome_domain_raw <- read.csv("../../data/Cb_bounds_df.csv", header = TRUE)

genome_domain <- genome_domain_raw %>%
  mutate(
    category = case_when(
      grepl("tip$", sub_region)    ~ "Tip",
      grepl("arm$", sub_region)    ~ "Arm",
      sub_region == "center"       ~ "Center",
      TRUE                         ~ NA_character_
    ),
    xmin = start / 1e6,
    xmax = stop  / 1e6
  ) %>%
  filter(!is.na(category)) %>%
  group_by(chrom) %>%
  mutate(
    xmin_norm = (xmin - min(xmin, na.rm = TRUE)) / (max(xmax, na.rm = TRUE) - min(xmin, na.rm = TRUE)),
    xmax_norm = (xmax - min(xmin, na.rm = TRUE)) / (max(xmax, na.rm = TRUE) - min(xmin, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = c("Tip", "Arm", "Center")))


# —————————————————————————————————————————
# （5）plot with domain shading
# —————————————————————————————————————————
make_stat_panel <- function(data, stat_filter) {
  ggplot() +
    # domain 
    geom_rect(
      data = genome_domain,
      inherit.aes = FALSE,
      aes(
        xmin = xmin_norm,
        xmax = xmax_norm,
        ymin = -Inf,
        ymax = +Inf,
        fill = category
      ),
      alpha = 0.4
    ) +
    geom_point(
      data = data %>% filter(stat_type == stat_filter),
      aes(x = x_norm, y = stat),
      color = "gray30", size = 0.2, alpha = 0.8, shape = 16
    ) +
    geom_smooth(
      data = data %>% filter(stat_type == stat_filter),
      aes(x = x_norm, y = stat),
      method = "loess", se = FALSE, span = 0.3, color = "lightgray"
    ) +
    facet_grid(region ~ chrom, scales = "fixed") +
    scale_fill_manual(values = genome_domain_colors) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Normalized genome position", y = stat_filter) +
    theme_bw() +
    theme(
      panel.grid        = element_blank(),
      legend.position   = "none",
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = 7, color = "black"),
      axis.title.y      = element_text(face = "bold", size = 8),
      axis.title.x      = element_text(face = "bold", size = 8),
      axis.text.x       = element_text(size = 6),
      axis.text.y       = element_text(size = 6)
    )
}





### add sample size 
Cb_indep_isotype_info_geo<-read.csv("../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv")

n_Australia<-Cb_indep_isotype_info_geo %>% 
  filter(geo == "Australia") %>% 
  nrow()
n_Australia
#57

isotype_by_lineage_from_Nic_raw<-read.table("../../data/From_Nic/isotype_byLineage_GeoLocAdmCol.tsv",
                                            header = TRUE,
                                            sep = '\t',
                                            comment.char  = "")

isotype_by_lineage<-isotype_by_lineage_from_Nic_raw %>% 
select(isotype,Sublineage,Lineage) %>% 
  left_join(Cb_indep_isotype_info_geo,by = "isotype") 


n_AD<-isotype_by_lineage %>% 
  filter(Lineage == "AD") %>%
  filter(geo == "Australia") %>% 
  nrow()
n_AD
#32

n_Tropical<-isotype_by_lineage %>% 
  filter(Lineage == "Tropical") %>% 
  filter(geo == "Australia")%>% 
  nrow()
n_Tropical
#25

n_All_Australia<-Cb_indep_isotype_info_geo %>% 
  filter(geo == "Australia") %>% 
  nrow()
n_All_Australia
#57


# intersect(n_AD$isotype,n_Tropical$isotype)
#character(0)

# all_diversity<-all_diversity %>% 
  
all_diversity<-all_diversity %>% 
  mutate(region = ifelse(region =="All Australia", paste0("All Australia (", n_All_Australia, ")"),
                         ifelse(region =="AD",paste0("AD (", n_AD, ")"),
                                ifelse(region =="Tropical",paste0("Tropical (", n_Tropical, ")"),region
                                       )
                                )
                         )
         ) %>% 
  mutate(region = factor(region,
                         levels = c(
                           unique(region[stringr::str_detect(region, "^All Australia")]),
                           unique(region[!stringr::str_detect(region, "^All Australia")])
                         )))



p_pi    <- make_stat_panel(all_diversity, "Nucleotide diversity (\u03C0)")
p_theta <- make_stat_panel(all_diversity, "Watterson's \u03B8")
p_d     <- make_stat_panel(all_diversity, "Tajima's D")

ggsave("Australia_pi_windowed_diversity.png",    p_pi,    width = 7, height = 4.5, units = "in", dpi = 300)
ggsave("Australia_theta_windowed_diversity.png", p_theta, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave("Australia_d_windowed_diversity.png",     p_d,     width = 7, height = 4.5, units = "in", dpi = 300)

