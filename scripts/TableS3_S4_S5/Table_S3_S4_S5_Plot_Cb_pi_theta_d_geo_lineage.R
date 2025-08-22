rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(purrr)



isotype_geo_raw<-read.csv("../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv")

n_all <- data.frame(
  Lineage  = "All",
  Frequency = nrow(isotype_geo_raw),
  stringsAsFactors = FALSE
)

n_geo<-isotype_geo_raw %>% 
  group_by(geo) %>% 
  summarise(Frequency = n()) %>% 
  bind_rows(
    tibble(
      geo = "Cb_non_cosmopolitan",
      Frequency = sum(.$Frequency) -
        filter(., geo == "cosmopolitan")$Frequency
    )
  )%>% 
  filter(geo %in% c("Asia","Australia","Caribbean", 
                    "Central America",
                    # "cosmopolitan", 
                    "Hawaii",
                    "Cb_non_cosmopolitan",
                    "Pacific","South America","Taiwan")) %>% 
  rename(Lineage=geo) %>% 
  mutate(Lineage = ifelse(Lineage == "Central America","Central_America",
                          ifelse(Lineage == "South America", "South_America",Lineage)))


  

sub_region_raw<-read.table("../../data/05.07.21_cb_subregion.bed")


sub_region_df <- sub_region_raw %>%
  filter(V1 != "CHROM") %>%
  mutate(
    V2 = as.numeric(V2),
    V3 = as.numeric(V3)
  ) %>%
  rename(
    chrom        = V1,
    region_start = V2,
    region_end   = V3,
    subregion    = V4
  )



regions_geo <- c("Cb_non_cosmopolitan","Asia","Australia",
                 "Caribbean", "Central_America", 
                 # "cosmopolitan",
                 "Hawaii","Pacific",
                 "South_America","Taiwan")

base_dir <- "../../processed_data/pi_theta_d_geo"


########################################
######## By center, arm, tip ########
########################################

all_geo_summary_combos <- purrr::map_dfr(regions_geo, function(region_name) {
  # 1) add mid 
  wd <- read_csv(
    file.path(base_dir, region_name, "chromosome_windows_diversity.csv"),
    col_types = cols()
  ) %>%
    mutate(mid = (window_start + window_stop) / 2)
  
  # 2) calculate each subregion
  sub_sum <- wd %>%
    inner_join(sub_region_df, by = "chrom") %>%
    filter(mid >= region_start, mid <= region_end) %>%
    group_by(Geo       = region_name,
             Subregion = subregion,
             Stat      = stat_type) %>%
    summarise(Mean_value = mean(stat, na.rm = TRUE), .groups = "drop")
  
  # 2.5) define combo regions：center, arm, tip
  combos <- list(
    center = "center",
    arm    = c("left_arm", "right_arm"),
    tip    = c("left_tip", "right_tip")
  )
  domain_sum <- purrr::map_dfr(names(combos), function(dom) {
    wd %>%
      inner_join(sub_region_df, by = "chrom") %>%
      filter(mid >= region_start,
             mid <= region_end,
             subregion %in% combos[[dom]]) %>%
      group_by(Geo  = region_name,
               Stat = stat_type) %>%
      summarise(Mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      mutate(Subregion = dom) %>%
      select(Geo, Subregion, Stat, Mean_value)
  })
  
  # 3) calculate whole genome-wide
  genome_sum <- wd %>%
    group_by(Geo  = region_name,
             Stat = stat_type) %>%
    summarise(Mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
    mutate(Subregion = "full") %>%
    select(Geo, Subregion, Stat, Mean_value)
  
  # merge
  dplyr::bind_rows(sub_sum, domain_sum, genome_sum)
})


# long table
print(all_geo_summary_combos)






###### All 719 isotypes ######
regions_all<-c("All")
base_dir_all<-"../../processed_data/pi_theta_d"
all_isotypes_summary_combos <- purrr::map_dfr(regions_all, function(region_name) {
  # 1) add mid 
  wd <- read_csv(
    file.path(base_dir_all, "chromosome_windows_diversity.csv"),
    col_types = cols()
  ) %>%
    mutate(mid = (window_start + window_stop) / 2)
  
  # 2) calculate each subregion
  sub_sum <- wd %>%
    inner_join(sub_region_df, by = "chrom") %>%
    filter(mid >= region_start, mid <= region_end) %>%
    group_by(Geo       = region_name,
             Subregion = subregion,
             Stat      = stat_type) %>%
    summarise(Mean_value = mean(stat, na.rm = TRUE), .groups = "drop")
  
  # 2.5) define combo regions：center, arm, tip
  combos <- list(
    center = "center",
    arm    = c("left_arm", "right_arm"),
    tip    = c("left_tip", "right_tip")
  )
  domain_sum <- purrr::map_dfr(names(combos), function(dom) {
    wd %>%
      inner_join(sub_region_df, by = "chrom") %>%
      filter(mid >= region_start,
             mid <= region_end,
             subregion %in% combos[[dom]]) %>%
      group_by(Geo  = region_name,
               Stat = stat_type) %>%
      summarise(Mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      mutate(Subregion = dom) %>%
      select(Geo, Subregion, Stat, Mean_value)
  })
  
  # 3) calculate whole genome-wide
  genome_sum <- wd %>%
    group_by(Geo  = region_name,
             Stat = stat_type) %>%
    summarise(Mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
    mutate(Subregion = "full") %>%
    select(Geo, Subregion, Stat, Mean_value)
  
  # merge
  dplyr::bind_rows(sub_sum, domain_sum, genome_sum)
})

all_isotypes_summary_combos





all_groups_combos<-rbind(all_isotypes_summary_combos,
                         all_geo_summary_combos)




all_groups_combos<-all_groups_combos %>%
  # mutate(Mean_value = round(Mean_value,4)) %>% 
  unique()



all_wide_combos <- all_groups_combos %>%
  pivot_wider(
    names_from  = Stat,
    values_from = Mean_value
  ) %>% 
  mutate(Geo = ifelse(Geo=="North_America","North America",
                                 ifelse(Geo=="South_America","South America",
                                        ifelse(Geo=="Cb_non_cosmopolitan","Non-cosmopolitan", 
                                               Geo)))) %>% 
  select(-mis)




###### pi. #######
pi_wide_combos<-all_wide_combos %>% 
  select(-theta,-d) %>% 
  pivot_wider(
    names_from  = Subregion,
    values_from = pi
  ) %>% 
  select(Geo,center,arm,tip,full) %>% 
  mutate(`fold change arm/center`= arm/tip) %>%
  mutate(across(where(is.numeric), ~ signif(., 4)))

head(pi_wide_combos)
write.table(pi_wide_combos,
            "geo_pi_armVScenter.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)



####### theta #####
theta_wide_combos<-all_wide_combos %>% 
  select(-pi,-d) %>% 
  pivot_wider(
    names_from  = Subregion,
    values_from = theta
  ) %>% 
  select(Geo,center,arm,tip,full) %>% 
  mutate(`fold change arm/center`= arm/tip) %>%
  mutate(across(where(is.numeric), ~ signif(., 4)))

head(theta_wide_combos)
write.table(theta_wide_combos,
            "geo_theta_armVScenter.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)





###### d ###
d_wide_combos<-all_wide_combos %>% 
  select(-pi,-theta) %>% 
  pivot_wider(
    names_from  = Subregion,
    values_from = d
  ) %>% 
  select(Geo,center,arm,tip,full) %>% 
  mutate(`fold change arm/center`= arm/tip) %>%
  mutate(across(where(is.numeric), ~ signif(., 4)))

head(d_wide_combos)
write.table(d_wide_combos,
            "geo_d_arm_center.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

















