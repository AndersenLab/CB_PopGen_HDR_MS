rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

sub_region_raw <- read.table("../../data/briggsae_genome_files/05.07.21_cb_subregion.bed")
sub_region_df <- sub_region_raw %>%
  dplyr::filter(V1 != "CHROM") %>%
  dplyr::mutate(
    V2 = as.numeric(V2),
    V3 = as.numeric(V3)
  ) %>%
  dplyr::rename(
    chrom        = V1,
    region_start = V2,
    region_end   = V3,
    subregion    = V4
  )

region_paths <- list(
  AD ="../../processed_data/diversity_and_divergence/pi_theta_d_rg/AD/",
  KD ="../../processed_data/diversity_and_divergence/pi_theta_d_rg/KD/",
  TD1 = "../../processed_data/diversity_and_divergence/pi_theta_d_rg/TD1",
  Temperate="../../processed_data/diversity_and_divergence/pi_theta_d_rg/Temperate/",
  TH="../../processed_data/diversity_and_divergence/pi_theta_d_rg/TH/",
  Tropical="../../processed_data/diversity_and_divergence/pi_theta_d_rg/Tropical/"
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
      dplyr::filter(chrom %in% chrom_filter, stat_type == "pi") %>%
      dplyr::mutate(mid = (window_start + window_stop) / 2) %>%
      dplyr::inner_join(sub_region_df, by = "chrom") %>%
      dplyr::filter(mid >= region_start, mid <= region_end) %>%
      dplyr::mutate(domain = case_when(
        subregion == "center" ~ "center",
        subregion %in% c("left_arm", "right_arm") ~ "arm",
        TRUE ~ NA_character_
      )) %>%
      dplyr::filter(!is.na(domain)) %>%
      dplyr::group_by(domain) %>%
      dplyr::summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(type = "pi", chrom = chrom_label, region = region)
  }
  
  autosomes <- process(df, c("I","II","III","IV","V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  bind_rows(autosomes, chrX)
}

all_pi_results <- imap_dfr(region_paths, div_calc)


wide_pi_results <- all_pi_results %>%
  tidyr::unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  tidyr::pivot_wider(
    id_cols = region,
    names_from = chrom_domain,
    values_from = mean_value
  ) %>% 
  dplyr::mutate(Autosomes_arm=round(Autosomes_arm,4)) %>% 
  dplyr::mutate(Autosomes_center=round(Autosomes_center,4)) %>% 
  dplyr::mutate(X_arm=round(X_arm,4)) %>% 
  dplyr::mutate(X_center=round(X_center,4)) %>% 
  dplyr::rename(`Autosomes arm`=Autosomes_arm,
         `Autosomes center`=Autosomes_center,
         `ChromX arm`=X_arm,
         `ChromX center`=X_center,
         Region=region) %>% 
  dplyr::mutate(Stat = "pi") 

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
      dplyr::filter(chrom %in% chrom_filter, stat_type == "theta") %>%
      dplyr::mutate(mid = (window_start + window_stop) / 2) %>%
      dplyr::inner_join(sub_region_df, by = "chrom") %>%
      dplyr::filter(mid >= region_start, mid <= region_end) %>%
      dplyr::mutate(domain = case_when(
        subregion == "center" ~ "center",
        subregion %in% c("left_arm", "right_arm") ~ "arm",
        TRUE ~ NA_character_
      )) %>%
      dplyr::filter(!is.na(domain)) %>%
      dplyr::group_by(domain) %>%
      dplyr::summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(type = "theta", chrom = chrom_label, region = region)
  }
  
  autosomes <- process(df, c("I","II","III","IV","V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  bind_rows(autosomes, chrX)
}

all_theta_results <- imap_dfr(region_paths, div_calc)


wide_theta_results <- all_theta_results %>%
  tidyr::unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  tidyr::pivot_wider(
    id_cols = region,
    names_from = chrom_domain,
    values_from = mean_value
  ) %>% 
  dplyr::mutate(Autosomes_arm=round(Autosomes_arm,4)) %>% 
  dplyr::mutate(Autosomes_center=round(Autosomes_center,4)) %>% 
  dplyr::mutate(X_arm=round(X_arm,4)) %>% 
  dplyr::mutate(X_center=round(X_center,4)) %>% 
  dplyr::rename(`Autosomes arm`=Autosomes_arm,
         `Autosomes center`=Autosomes_center,
         `ChromX arm`=X_arm,
         `ChromX center`=X_center,
         Region=region) %>% 
  dplyr::mutate(Stat = "theta")


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
      dplyr::filter(chrom %in% chrom_filter, stat_type == "d") %>%
      dplyr::mutate(mid = (window_start + window_stop) / 2) %>%
      dplyr::inner_join(sub_region_df, by = "chrom") %>%
      dplyr::filter(mid >= region_start, mid <= region_end) %>%
      dplyr::mutate(domain = case_when(
        subregion == "center" ~ "center",
        subregion %in% c("left_arm", "right_arm") ~ "arm",
        TRUE ~ NA_character_
      )) %>%
      dplyr::filter(!is.na(domain)) %>%
      dplyr::group_by(domain) %>%
      dplyr::summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(type = "Tajima's D", chrom = chrom_label, region = region)
  }
  
  autosomes <- process(df, c("I","II","III","IV","V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  bind_rows(autosomes, chrX)
}

all_d_results <- imap_dfr(region_paths, div_calc)


wide_d_results <- all_d_results %>%
  tidyr::unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  tidyr::pivot_wider(
    id_cols = region,
    names_from = chrom_domain,
    values_from = mean_value
  ) %>% 
  dplyr::mutate(Autosomes_arm=round(Autosomes_arm,4)) %>% 
  dplyr::mutate(Autosomes_center=round(Autosomes_center,4)) %>% 
  dplyr::mutate(X_arm=round(X_arm,4)) %>% 
  dplyr::mutate(X_center=round(X_center,4)) %>% 
  dplyr::rename(`Autosomes arm`=Autosomes_arm,
         `Autosomes center`=Autosomes_center,
         `ChromX arm`=X_arm,
         `ChromX center`=X_center,
         Region=region) %>% 
  dplyr::mutate(Stat = "Tajima's D") 


merged_wide_table<-rbind(wide_pi_results,
                         wide_theta_results,
                         wide_d_results)


merged_wide_table_foldchange<-merged_wide_table %>% 
  dplyr::mutate(`fold change Autosomes arm/center` = (round(`Autosomes arm`/`Autosomes center`,2))) %>% 
  dplyr::mutate(`fold change ChromX arm/center` = (round(`ChromX arm`/`ChromX center`,2))) %>% 
  dplyr::relocate(Stat, .after = last_col())

write.table(merged_wide_table_foldchange,
            "../../supplementary_data/SD12_All_groups_merged_pi_theta_d_Autosomal_Xarm_arm_center_foldchange.tsv",
            quote = FALSE,
            row.names = FALSE,
            sep = '\t')

