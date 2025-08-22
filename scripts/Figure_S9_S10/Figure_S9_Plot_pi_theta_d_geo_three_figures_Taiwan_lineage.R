rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


source("../utilities.R")

# （1）all isotypes of the geo 
geo_df <- read.csv(
  "../../processed_data/pi_theta_d_geo/Taiwan/chromosome_windows_diversity.csv",
  header = TRUE
) %>%
  filter(stat_type %in% c("pi", "theta", "d")) %>%
  mutate(
    region    = "All Taiwan",
    stat_type = recode(stat_type,
                       pi    = "Nucleotide diversity (\u03C0)",
                       theta = "Watterson's \u03B8",
                       d     = "Tajima's D")
  ) %>%
  group_by(chrom) %>%
  mutate(x= x/1e6) %>% 
  # mutate(
  #   x_norm = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  # ) %>%
  ungroup()


# （2）lineage function
read_window_data <- function(region_name,
                             base_dir = "../../processed_data/pi_theta_d_by_lineage/Taiwan_lineage") {
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
    mutate(x= x/1e6) %>% 
    # mutate(
    #   x_norm = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    # ) %>%
    ungroup()
}

# （3）read lineages
regions <- c("TD1", "TH", "Tropical")
all_diversity <- lapply(regions, read_window_data) %>% bind_rows()

# merge
all_diversity <- bind_rows(geo_df, all_diversity) %>%
  mutate(
    region = factor(region,
                    levels = c("All Taiwan", "TD1", "TH", "Tropical")),
    stat_type = factor(stat_type,
                       levels = c("Nucleotide diversity (π)",
                                  "Watterson's θ",
                                  "Tajima's D"))
  )

# all_diversity <- bind_rows(geo_df,all_diversity) %>%
#   mutate(
#     region = factor(region,
#                     levels = c("All Taiwan", regions)),
#     stat_type = factor(stat_type,
#                        levels = c("Nucleotide diversity (\u03C0)",
#                                   "Watterson's \u03B8",
#                                   "Tajima's D"))
#   )




### (4) read genome domain
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


# (5) plot
make_stat_panel <- function(data, stat_filter) {
  ggplot() +
    
    geom_rect(
      data = genome_domain,
      inherit.aes = FALSE,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = +Inf,
        fill = category
      ),
      alpha = 0.5
    ) +
    
    geom_point(
      data = data %>% filter(stat_type == stat_filter),
      aes(x = x, y = stat),
      color = "gray30", size = 0.2, alpha = 0.8, shape = 16
    ) +
    geom_smooth(
      data = data %>% filter(stat_type == stat_filter),
      aes(x = x, y = stat),
      method = "loess", se = FALSE, span = 0.3, color = "lightgray"
    ) +
    
    facet_grid(region ~ chrom, scales = "free_x") +
    
    scale_fill_manual(values = genome_domain_colors) +
    labs(x = "Physical genome position (Mb)",
         y = stat_filter) +
    theme_bw() +
    theme(
      panel.grid       = element_blank(),
      legend.position  = "none",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 7),
      axis.title       = element_text(face = "bold", size = 8),
      axis.text        = element_text(size = 6)
    )
}







### add sample size 
Cb_indep_isotype_info_geo<-read.csv("../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv")

isotype_by_lineage_from_Nic_raw<-read.table("../../data/From_Nic/isotype_byLineage_GeoLocAdmCol.tsv",
                                            header = TRUE,
                                            sep = '\t',
                                            comment.char  = "")

isotype_by_lineage<-isotype_by_lineage_from_Nic_raw %>% 
  select(isotype,Sublineage,Lineage) %>% 
  left_join(Cb_indep_isotype_info_geo,by = "isotype") 


n_TD1<-isotype_by_lineage %>% 
  filter(Lineage == "TD1") %>%
  filter(geo == "Taiwan") %>% 
  nrow()
n_TD1
#27


n_TH<-isotype_by_lineage %>% 
  filter(Lineage == "TH") %>%
  filter(geo == "Taiwan") %>% 
  nrow()
n_TH
#76


n_Tropical<-isotype_by_lineage %>% 
  filter(Lineage == "Tropical") %>% 
  filter(geo == "Taiwan")%>% 
  nrow()
n_Tropical
#41


n_All_Taiwan<-Cb_indep_isotype_info_geo %>% 
  filter(geo == "Taiwan") %>% 
  nrow()
n_All_Taiwan
#153


# intersect(n_AD$isotype,n_Tropical$isotype)
#character(0)

# all_diversity<-all_diversity %>% 

all_diversity<-all_diversity %>% 
  mutate(region = ifelse(region =="All Taiwan", paste0("All Taiwan (", n_All_Taiwan, ")"),
                         ifelse(region =="TD1",paste0("TD1 (", n_TD1, ")"),
                                ifelse(region =="TH",paste0("TH (", n_TH, ")"),
                                       ifelse(region =="Tropical",paste0("Tropical (", n_Tropical, ")"),region
                                              )
                                )
                         )
  )
  ) %>% 
  mutate(region = factor(region,
                         levels = c(
                           unique(region[stringr::str_detect(region, "^All Taiwan")]),
                           unique(region[!stringr::str_detect(region, "^All Taiwan")])
                         )))









p_pi    <- make_stat_panel(all_diversity, "Nucleotide diversity (\u03C0)")
p_theta <- make_stat_panel(all_diversity, "Watterson's \u03B8")
p_d     <- make_stat_panel(all_diversity, "Tajima's D")


ggsave(
  filename = "Taiwan_pi_windowed_diversity.png",
  plot     = p_pi,
  device   = "png",
  width    = 7.5,
  height   = 6,
  units    = "in",
  dpi      = 300
)



ggsave(
  filename = "Taiwan_theta_windowed_diversity.png",
  plot     = p_theta,
  device   = "png",
  width    = 7.5,
  height   = 6,
  units    = "in",
  dpi      = 300
)



ggsave(
  filename = "Taiwan_d_windowed_diversity.png",
  plot     = p_d,
  device   = "png",
  width    = 7.5,
  height   = 6,
  units    = "in",
  dpi      = 300
)

