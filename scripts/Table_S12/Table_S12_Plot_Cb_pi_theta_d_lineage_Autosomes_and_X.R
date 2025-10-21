rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

sub_region_raw <- read.table("../../data/05.07.21_cb_subregion.bed")
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

region_paths <- list(
  AD ="../../processed_data/pi_theta_d_by_lineage/Australia_lineage/AD/",
  KD ="../../processed_data/pi_theta_d_TH_KD_TD2/KD_lineage/KD/",
  TD1 = "../../processed_data/pi_theta_d_by_lineage/Taiwan_lineage/TD1",
  Temperate="../../processed_data/Temperate_hard_filtered_and_LD_pruned/Cb_hard_filtered_Temperate/",
  TH="../../processed_data/TH_hard_filtered_pi_theta_d/Cb_hard_filtered_TH/",
  Tropical="../../processed_data/Tropical_hard_filtered_and_LD_pruned/Cb_hard_filtered_Tropical/"
)

div_calc <- function(path,region) {
  file_path <- file.path(path, "chromosome_windows_diversity.csv")
  
  message("Looking for file: ", normalizePath(file_path, mustWork = FALSE))
  
  if (!file.exists(file_path)) {
    warning(paste("Missing file:", file_path))
    return(NULL)
  }
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  process <- function(df, chrom_filter, chrom_label) {
    df %>%
      filter(chrom %in% chrom_filter, stat_type == "pi") %>%
      mutate(mid = (window_start + window_stop) / 2) %>%
      inner_join(sub_region_df, by = "chrom") %>%
      filter(mid >= region_start, mid <= region_end) %>%
      mutate(domain = case_when(
        subregion == "center" ~ "center",
        subregion %in% c("left_arm", "right_arm") ~ "arm",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(domain)) %>%
      group_by(domain) %>%
      summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      mutate(type = "pi", chrom = chrom_label, region = region)
  }
  
  autosomes <- process(df, c("I","II","III","IV","V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  bind_rows(autosomes, chrX)
}

all_pi_results <- imap_dfr(region_paths, div_calc)


wide_pi_results <- all_pi_results %>%
  unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  pivot_wider(
    id_cols = region,
    names_from = chrom_domain,
    values_from = mean_value
  ) %>% 
  mutate(Autosomes_arm=round(Autosomes_arm,4)) %>% 
  mutate(Autosomes_center=round(Autosomes_center,4)) %>% 
  mutate(X_arm=round(X_arm,4)) %>% 
  mutate(X_center=round(X_center,4)) %>% 
  rename(`Autosomes arm`=Autosomes_arm,
         `Autosomes center`=Autosomes_center,
         `ChromX arm`=X_arm,
         `ChromX center`=X_center,
         Region=region) %>% 
  mutate(Stat = "pi") 

div_calc <- function(path,region) {
  file_path <- file.path(path, "chromosome_windows_diversity.csv")
  
  message("Looking for file: ", normalizePath(file_path, mustWork = FALSE))
  
  if (!file.exists(file_path)) {
    warning(paste("Missing file:", file_path))
    return(NULL)
  }
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  process <- function(df, chrom_filter, chrom_label) {
    df %>%
      filter(chrom %in% chrom_filter, stat_type == "theta") %>%
      mutate(mid = (window_start + window_stop) / 2) %>%
      inner_join(sub_region_df, by = "chrom") %>%
      filter(mid >= region_start, mid <= region_end) %>%
      mutate(domain = case_when(
        subregion == "center" ~ "center",
        subregion %in% c("left_arm", "right_arm") ~ "arm",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(domain)) %>%
      group_by(domain) %>%
      summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      mutate(type = "theta", chrom = chrom_label, region = region)
  }
  
  autosomes <- process(df, c("I","II","III","IV","V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  bind_rows(autosomes, chrX)
}

all_theta_results <- imap_dfr(region_paths, div_calc)


wide_theta_results <- all_theta_results %>%
  unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  pivot_wider(
    id_cols = region,
    names_from = chrom_domain,
    values_from = mean_value
  ) %>% 
  mutate(Autosomes_arm=round(Autosomes_arm,4)) %>% 
  mutate(Autosomes_center=round(Autosomes_center,4)) %>% 
  mutate(X_arm=round(X_arm,4)) %>% 
  mutate(X_center=round(X_center,4)) %>% 
  rename(`Autosomes arm`=Autosomes_arm,
         `Autosomes center`=Autosomes_center,
         `ChromX arm`=X_arm,
         `ChromX center`=X_center,
         Region=region) %>% 
  mutate(Stat = "theta")


div_calc <- function(path,region) {
  file_path <- file.path(path, "chromosome_windows_diversity.csv")
  
  message("Looking for file: ", normalizePath(file_path, mustWork = FALSE))
  
  if (!file.exists(file_path)) {
    warning(paste("Missing file:", file_path))
    return(NULL)
  }
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  process <- function(df, chrom_filter, chrom_label) {
    df %>%
      filter(chrom %in% chrom_filter, stat_type == "d") %>%
      mutate(mid = (window_start + window_stop) / 2) %>%
      inner_join(sub_region_df, by = "chrom") %>%
      filter(mid >= region_start, mid <= region_end) %>%
      mutate(domain = case_when(
        subregion == "center" ~ "center",
        subregion %in% c("left_arm", "right_arm") ~ "arm",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(domain)) %>%
      group_by(domain) %>%
      summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      mutate(type = "Tajima's D", chrom = chrom_label, region = region)
  }
  
  autosomes <- process(df, c("I","II","III","IV","V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  bind_rows(autosomes, chrX)
}

all_d_results <- imap_dfr(region_paths, div_calc)


wide_d_results <- all_d_results %>%
  unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  pivot_wider(
    id_cols = region,
    names_from = chrom_domain,
    values_from = mean_value
  ) %>% 
  mutate(Autosomes_arm=round(Autosomes_arm,4)) %>% 
  mutate(Autosomes_center=round(Autosomes_center,4)) %>% 
  mutate(X_arm=round(X_arm,4)) %>% 
  mutate(X_center=round(X_center,4)) %>% 
  rename(`Autosomes arm`=Autosomes_arm,
         `Autosomes center`=Autosomes_center,
         `ChromX arm`=X_arm,
         `ChromX center`=X_center,
         Region=region) %>% 
  mutate(Stat = "Tajima's D") 










merged_wide_table<-rbind(wide_pi_results,
                         wide_theta_results,
                         wide_d_results)


write.table(merged_wide_table,
            "All_groups_merged_pi_theta_d_Autosomal_Xarm_arm_center.tsv",
            quote = FALSE,
            row.names = FALSE,
            sep = '\t')








merged_wide_table_foldchange<-merged_wide_table %>% 
  mutate(`fold change Autosomes arm/center` = (round(`Autosomes arm`/`Autosomes center`,2))) %>% 
  mutate(`fold change ChromX arm/center` = (round(`ChromX arm`/`ChromX center`,2))) %>% 
  relocate(Stat, .after = last_col())

write.table(merged_wide_table_foldchange,
            "All_groups_merged_pi_theta_d_Autosomal_Xarm_arm_center_foldchange.tsv",
            quote = FALSE,
            row.names = FALSE,
            sep = '\t')

