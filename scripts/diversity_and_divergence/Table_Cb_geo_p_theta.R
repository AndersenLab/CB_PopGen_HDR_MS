

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)


calculate_diversity<-function(custom_region,
                              diversity_type,
                              custom_path = NULL){

  if (!is.null(custom_path)) {
    Geo_stat_raw <- read.csv(custom_path)
  } else {
  file_path <- c('../../processed_data/diversity_and_divergence/pi_theta_d_geo/region/chromosome_windows_diversity.csv')
  new_file_path <- gsub("region", custom_region, file_path)
  Geo_stat_raw <- read.csv(new_file_path)
  }


diversity <- Geo_stat_raw %>%
  dplyr::filter(stat_type == diversity_type) %>%
  dplyr::select(stat) %>% 
  na.omit()
  
output<-mean(diversity$stat)
return(output)
    

}

Cb_cosmopolitan_isotype_removed_theta<-calculate_diversity("Cb_non_cosmopolitan","theta")
Africa_theta<-calculate_diversity("Africa","theta")
Asia_theta<-calculate_diversity("Asia","theta")
Australia_theta<-calculate_diversity("Australia","theta")
Caribbean_theta<-calculate_diversity("Caribbean","theta")
Central_America_theta<-calculate_diversity("Central_America","theta")
Europe_theta<-calculate_diversity("Europe","theta")
Hawaii_theta<-calculate_diversity("Hawaii","theta")
North_America_theta<-calculate_diversity("North_America","theta")
Pacific_theta<-calculate_diversity("Pacific","theta")
South_America_theta<-calculate_diversity("South_America","theta")
Taiwan_theta<-calculate_diversity("Taiwan","theta")
Ce_all_theta<-calculate_diversity(diversity_type="theta",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_pi_theta_d/chromosome_windows_diversity.csv")
Ce_Hawaii_theta<-calculate_diversity(diversity_type="theta",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_Hawaii/chromosome_windows_diversity.csv")
Cb_Thomas_theta<-calculate_diversity(diversity_type="theta",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Thomas_pi_theta_d/chromosome_windows_diversity.csv")


Cb_cosmopolitan_isotype_removed_pi<-calculate_diversity("Cb_non_cosmopolitan","pi")
Africa_pi<-calculate_diversity("Africa","pi")
Asia_pi<-calculate_diversity("Asia","pi")
Australia_pi<-calculate_diversity("Australia","pi")
Caribbean_pi<-calculate_diversity("Caribbean","pi")
Central_America_pi<-calculate_diversity("Central_America","pi")
Europe_pi<-calculate_diversity("Europe","pi")
Hawaii_pi<-calculate_diversity("Hawaii","pi")
North_America_pi<-calculate_diversity("North_America","pi")
Pacific_pi<-calculate_diversity("Pacific","pi")
South_America_pi<-calculate_diversity("South_America","pi")
Taiwan_pi<-calculate_diversity("Taiwan","pi")
Ce_all_pi<-calculate_diversity(diversity_type="pi",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_pi_theta_d/chromosome_windows_diversity.csv")
Ce_Hawaii_pi<-calculate_diversity(diversity_type="pi",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_Hawaii/chromosome_windows_diversity.csv")
Cb_Thomas_pi<-calculate_diversity(diversity_type="pi",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Thomas_pi_theta_d/chromosome_windows_diversity.csv")


Cb_cosmopolitan_isotype_removed_d<-calculate_diversity("Cb_non_cosmopolitan","d")
Africa_d<-calculate_diversity("Africa","d")
Asia_d<-calculate_diversity("Asia","d")
Australia_d<-calculate_diversity("Australia","d")
Caribbean_d<-calculate_diversity("Caribbean","d")
Central_America_d<-calculate_diversity("Central_America","d")
Europe_d<-calculate_diversity("Europe","d")
Hawaii_d<-calculate_diversity("Hawaii","d")
North_America_d<-calculate_diversity("North_America","d")
Pacific_d<-calculate_diversity("Pacific","d")
South_America_d<-calculate_diversity("South_America","d")
Taiwan_d<-calculate_diversity("Taiwan","d")
Ce_all_d<-calculate_diversity(diversity_type="d",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_pi_theta_d/chromosome_windows_diversity.csv")
Ce_Hawaii_d<-calculate_diversity(diversity_type="d",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_Hawaii/chromosome_windows_diversity.csv")
Cb_Thomas_d<-calculate_diversity(diversity_type="d",custom_path = "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Thomas_pi_theta_d/chromosome_windows_diversity.csv")



species_wide_stats<-read.csv(file = "../../processed_data/diversity_and_divergence/pi_theta_d/chromosome_windows_diversity.csv")

species_pi <- species_wide_stats %>%
  dplyr::filter(stat_type == "pi") %>%
  dplyr::select(stat) %>% 
  na.omit()
species_pi<-mean(species_pi$stat)

species_theta <- species_wide_stats %>%
  dplyr::filter(stat_type == "theta") %>%
  dplyr::select(stat) %>% 
  na.omit()
species_theta<-mean(species_theta$stat)

species_d <- species_wide_stats %>%
  dplyr::filter(stat_type == "d") %>%
  dplyr::select(stat) %>% 
  na.omit()
species_d<-mean(species_d$stat)


table_geo_p_theta_d <- data.frame(
  region = rep(c("All",
                 "Non-cosmopolitan",
                 "Asia", "Australia", "Caribbean", "Central_America", 
                 "Hawaii", 
                 "Pacific","South_America", "Taiwan",
                 "C. elegans - All", "C. elegans - Hawaii",
                 "C. briggsae - Thomas"), each = 3),
  stat = rep(c("theta", "pi","d"), times = 13),
  values = c(species_theta, species_pi, species_d, 
             Cb_cosmopolitan_isotype_removed_theta,Cb_cosmopolitan_isotype_removed_pi, Cb_cosmopolitan_isotype_removed_d,
             Asia_theta, Asia_pi, Asia_d,
             Australia_theta, Australia_pi, Australia_d, 
             Caribbean_theta, Caribbean_pi, Caribbean_d,
             Central_America_theta, Central_America_pi, Central_America_d,
             Hawaii_theta, Hawaii_pi, Hawaii_d,
             Pacific_theta, Pacific_pi, Pacific_d,
             South_America_theta, South_America_pi, South_America_d,
             Taiwan_theta, Taiwan_pi, Taiwan_d,
             Ce_all_theta,Ce_all_pi,Ce_all_d,
             Ce_Hawaii_theta, Ce_Hawaii_pi, Ce_Hawaii_d,
             Cb_Thomas_theta, Cb_Thomas_pi, Cb_Thomas_d)
)


table_geo_p_theta_d$values<-round(table_geo_p_theta_d$values,4)

table_geo_p_theta_d<-table_geo_p_theta_d %>% 
  dplyr::arrange(stat)



geo_diversity_result <- table_geo_p_theta_d %>%
  tidyr::pivot_wider(names_from = stat, values_from = values) %>%
  dplyr::filter(!is.na(pi) & !is.na(theta)& !is.na(d)) %>%
  dplyr::select(region, pi, theta, d) %>% 
  dplyr::rename(Region=region,"\u03C0 value"=pi,"\u03B8 value"=theta, "Tajima's D"=d)

geo_freq <- read.csv(file = "../../processed_data/Geo_info/Cb_isotype_geo_freq.csv")

geo_freq_col<-geo_freq %>%
  dplyr::mutate(geo = ifelse(geo == "South America", "South_America", geo)) %>% 
  dplyr::mutate(geo = ifelse(geo == "Central America", "Central_America", geo)) %>% 
  dplyr::mutate(geo = ifelse(geo == "North America", "North_America", geo)) 

geo_freq_col<-geo_freq_col %>% 
  filter(!(geo %in% "unknown"))
geo_freq_col<-rbind(geo_freq_col, c("Non-cosmopolitan", 715-(geo_freq_col %>%
                                      filter(geo == "Cosmopolitan") %>% 
                                      pull(frequency)))) 
geo_freq_col<-rbind(geo_freq_col, c("All", 715)) 


geo_freq_col$geo <- factor(geo_freq_col$geo, 
                           levels = c("All",
                                      "Non-cosmopolitan","Africa", 
                                      "Asia", "Australia", "Caribbean", "Central_America", 
                                      "Europe", "Hawaii", "North_America",
                                      "Pacific","South_America", "Taiwan"))
geo_freq_col_tmp1<-geo_freq_col %>% 
  dplyr::arrange(geo) %>% 
  dplyr::rename(Region=geo,"Number of strains"=frequency) %>% 
  mutate(`Number of strains` = as.integer(`Number of strains`)) 




Ce_isotype_groups_raw<-readr::read_tsv("../../data/Ce/isotype_groups.tsv")
Ce_isotypes_all<-Ce_isotype_groups_raw %>% 
  select(strain,isotype) %>% 
  filter(strain == isotype)
nrow(Ce_isotypes_all)
Ce_all_n<-nrow(Ce_isotypes_all)

Ce_Hawaii_list<-read.csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Ce_Hawaii/Ce_Hawaii_isotypes.csv")
nrow(Ce_Hawaii_list)
Ce_Hawaii_n<-nrow(Ce_Hawaii_list)

Cb_Thomas_n<-36



geo_freq_col_final <- geo_freq_col_tmp1 %>%
  add_row(
    Region = "C. elegans - All",
    `Number of strains` = Ce_all_n
  ) %>% 
  add_row(
    Region = "C. elegans - Hawaii",
    `Number of strains` = Ce_Hawaii_n
  ) %>% 
  add_row(
    Region = "C. briggsae - Thomas",
    `Number of strains` = Cb_Thomas_n) 


all_result<-dplyr::left_join(geo_diversity_result,
                             geo_freq_col_final,
                             by = c("Region"))


# write.csv(file = "Table_geo_p_theta_d.csv",all_result,quote = FALSE,
#           row.names = FALSE)

write.csv(file = "../../supplementary_data/SD4_geo_p_theta.csv",all_result,quote = FALSE,
          row.names = FALSE)








