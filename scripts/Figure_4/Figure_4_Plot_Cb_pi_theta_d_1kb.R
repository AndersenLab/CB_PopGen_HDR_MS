rm(list = ls())


library(dplyr)
library(ggplot2)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(8,"Set3")
brewer.pal(8,"Set3")[4:6]
display.brewer.pal(9,"Greys")
brewer.pal(9,"Greys")


source("../utilities.R")


species_wide_diversity_raw<-read.csv("../../processed_data/pi_theta_d_1kb/sub_region_diversity.csv",
         header = TRUE,
         row.names = 1)

species_wide_diversity<- species_wide_diversity_raw %>% 
  dplyr::mutate(sub_region = ifelse(sub_region == "left_arm", "left arm", sub_region)) %>% 
  dplyr::mutate(sub_region = ifelse(sub_region == "left_tip", "left tip", sub_region)) %>% 
  dplyr::mutate(sub_region = ifelse(sub_region == "right_arm", "right arm", sub_region)) %>% 
  dplyr::mutate(sub_region = ifelse(sub_region == "right_tip", "right tip", sub_region))
  
species_wide_diversity$sub_region <- factor(species_wide_diversity$sub_region, 
                                            levels = c("left tip", "left arm", "center", "right arm", "right tip", "full"))



colnames(species_wide_diversity)

# pi 
plot_pi<-species_wide_diversity %>%
  dplyr::select(chrom,sub_region,pi) 

# theta
plot_theta<-species_wide_diversity %>%
  dplyr::select(chrom,sub_region,theta)
# d
plot_d<-species_wide_diversity %>%
  dplyr::select(chrom,sub_region,d)




# plot pi heatmap
p_pi <- ggplot(plot_pi, aes(x = sub_region, y = chrom, fill = pi)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient(low = "lightyellow", high = "orange") +  
  ggplot2::theme_minimal() +  
  # ggplot2::labs(x = "Region", y = "Population") +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
  # theme(legend.position = c(0.2, 0.7))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(t = -0.1, unit = "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6, angle = 45, hjust = 1, vjust = 1.4))+
  ggplot2::ggtitle(expression("Nucleotide diversity (" * pi * ")"))+
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 8, vjust = -2, color = "black"),
                 axis.text.x = element_text(size = 6, color = "black"),  
                 axis.text.y = element_text(size = 6, color = "black"))+
  # ggplot2::ggtitle("Nucleotide diversity (\u03C0)")+
  ggplot2::labs(x = NULL, y = NULL)+
  ggplot2::geom_text(aes(label = round(pi, 5)), color = "black", size = 1.9)+  # add text labels
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
p_pi






# plot theta heatmap
p_theta <- ggplot(plot_theta, aes(x = sub_region, y = chrom, fill = theta)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient(low = "lightyellow", high = "orange") +  
  ggplot2::theme_minimal() +  
  # ggplot2::labs(x = "Region", y = "Population") +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
  # theme(legend.position = c(0.2, 0.7))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(t = -0.1, unit = "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6, angle = 45, hjust = 1, vjust = 1.4))+
  ggplot2::ggtitle(expression("Watterson's " * theta * ""))+
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 8, vjust = -2, color = "black"),
                 axis.text.x = element_text(size = 6, color = "black"),  
                 axis.text.y = element_text(size = 6, color = "black"))+
  # ggplot2::ggtitle("Nucleotide diversity (\u03C0)")+
  ggplot2::labs(x = NULL, y = NULL)+
  ggplot2::geom_text(aes(label = round(theta, 5)), color = "black", size = 1.9)+  # add text labels
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
p_theta







# plot d heatmap
p_d <- ggplot(plot_d, aes(x = sub_region, y = chrom, fill = d)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient(low = "lightyellow", high = "orange") +  
  ggplot2::theme_minimal() +  
  # ggplot2::labs(x = "Region", y = "Population") +
  ggplot2::theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
  # theme(legend.position = c(0.2, 0.7))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(t = -0.1, unit = "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6, angle = 45, hjust = 1, vjust = 1.4))+
  ggplot2::ggtitle(expression("Tajima's D"))+
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 8, vjust = -2, color = "black"),
                 axis.text.x = element_text(size = 6, color = "black"),  
                 axis.text.y = element_text(size = 6, color = "black"))+
  # ggplot2::ggtitle("Nucleotide diversity (\u03C0)")+
  ggplot2::labs(x = NULL, y = NULL)+
  ggplot2::geom_text(aes(label = round(d, 3)), color = "black", size = 1.9)+  # add text labels
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )
p_d


# plot a-c
library(ggpubr)
p_a_c <- ggpubr::ggarrange(p_pi, p_theta, p_d, 
                           ncol = 3,       
                           nrow = 1,        
                           labels = c("a", "b", "c"))     
p_a_c

ggsave("Cb_p_a_c.pdf", plot = p_a_c, width = 7.5, height = 2.5, units = "in")







###################################
######### genome wide loess #######
###################################


# Input data was kindly generated by Ryan


# install.packages("extrafont")
# font_import("Arial")

# rm(list = ls())

library(extrafont)
library(dplyr)
library(ggplot2)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(8,"Set3")
brewer.pal(8,"Set3")[4:6]
display.brewer.pal(9,"Greys")
brewer.pal(9,"Greys")

window_diversity<-read.csv("../../processed_data/pi_theta_d_1kb/chromosome_windows_diversity.csv",
                           header = TRUE,
                           row.names = 1)


# calculate mean theta and mean pi
window_diversity_pi<-window_diversity%>%
  dplyr::filter(stat_type == "pi") %>% 
  na.omit()
mean_pi <- mean(window_diversity_pi$stat)
mean_pi

window_diversity_theta<-window_diversity%>%
  dplyr::filter(stat_type == "theta") %>% 
  na.omit()
mean_theta <- mean(window_diversity_theta$stat)
mean_theta




###  read genome domain
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




#### Plot all windowed diversity stats
windowed_div_stats <- function(windows_df){
  diversity <- windows_df %>%
    filter(stat_type %in% c("pi","theta","d")) %>%
    mutate(
      stat_type = case_when(
        stat_type == "pi"    ~ "Nucleotide diversity (\u03C0)",
        stat_type == "theta" ~ "Watterson's \u03B8",
        stat_type == "d"     ~ "Tajima's D"
      ),
      x_mb = x / 1e6
    )
  
  diversity$stat_type <- factor(
    diversity$stat_type,
    levels = c("Nucleotide diversity (\u03C0)", 
               "Watterson's \u03B8", 
               "Tajima's D")
  )
  
  ggplot() +
    geom_rect(
      data = genome_domain,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    geom_point(
      data = diversity,
      aes(x = x_mb, y = stat),
      color = "gray30", size = 0.05, alpha = 0.8, shape = 16
    ) +
    geom_smooth(
      data = diversity,
      aes(x = x_mb, y = stat),
      method = "loess", se = FALSE, span = 0.3, color = "lightgray"
    ) +
    facet_grid(stat_type ~ chrom, scales = "free") +
    scale_fill_manual(values = genome_domain_colors) +
    xlab("Physical genome position (Mb)") +
    ylab("Diversity statistic") +
    theme_bw() +
    theme(
      legend.position      = "none",
      panel.grid           = element_blank(),
      strip.background     = element_blank(),
      strip.text           = element_text(face = "bold", size = 7, color = "#525252"),
      axis.title           = element_text(face = "bold", size = 9),
      axis.text            = element_text(size = 6)
    )
}




result_plots <- windowed_div_stats(window_diversity)
result_plots

# save image 4.5 * 7 inches
ggsave("raw_Cb_pi_theta_d_2025.pdf", plot = result_plots, width = 7, height = 4.5, units = "in")








### merge all
library(ggpubr)

plot_all<-ggarrange(
  ggpubr::ggarrange(p_pi, p_theta, p_d, nrow = 1, labels = c("a", "b", "c")),
  result_plots,
  heights = c(0.357, 0.643),  
  labels = c("","d"),  
  nrow = 2  
)
plot_all

ggsave("raw_Plot_pi_theta_d.pdf", plot = plot_all, width = 7, height = 7, units = "in")
ggsave("Plot_pi_theta_d.png", plot = plot_all, width = 7, height = 7, units = "in", dpi = 600)
















#################################################
####### Pi - Autosomes vs. Chromosome X ######
sub_region_raw<-read.table("../../data/05.07.21_cb_subregion.bed")

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

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



### Autosomes ##

Autosomes_pi_all<-window_diversity %>% 
  filter(chrom %in% c("I","II","III","IV","V")) %>% 
  filter(stat_type =="pi") %>% 
  mutate(mid = (window_start + window_stop) / 2) 



Autosomes_pi_with_region <- Autosomes_pi_all %>%
  inner_join(sub_region_df, by = "chrom") %>%
  filter(mid >= region_start, mid <= region_end)

Autosomes_pi_with_region <- Autosomes_pi_with_region %>%
  mutate(domain = case_when(
    subregion == "center"                 ~ "center",
    subregion %in% c("left_arm","right_arm") ~ "arm",
    TRUE                                  ~ NA_character_
  )) %>%
  filter(!is.na(domain))

Autosomes_mean_pi_by_chr <- Autosomes_pi_with_region %>%
  group_by(domain) %>%
  summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>% 
  mutate(type ="pi") %>% 
  mutate(chrom = "Autosomes")


Autosomes_fold_pi<-Autosomes_mean_pi_by_chr$mean_value[1]/Autosomes_mean_pi_by_chr$mean_value[2]
Autosomes_fold_pi
### 1.600007






### Chromosome_X ##

Chromosome_X_pi_all<-window_diversity %>% 
  filter(chrom %in% c("X")) %>% 
  filter(stat_type =="pi") %>% 
  mutate(mid = (window_start + window_stop) / 2) 



Chromosome_X_pi_with_region <- Chromosome_X_pi_all %>%
  inner_join(sub_region_df, by = "chrom") %>%
  filter(mid >= region_start, mid <= region_end)

Chromosome_X_pi_with_region <- Chromosome_X_pi_with_region %>%
  mutate(domain = case_when(
    subregion == "center"                 ~ "center",
    subregion %in% c("left_arm","right_arm") ~ "arm",
    TRUE                                  ~ NA_character_
  )) %>%
  filter(!is.na(domain))

Chromosome_X_mean_pi_by_chr <- Chromosome_X_pi_with_region %>%
  group_by(domain) %>%
  summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>% 
  mutate(type ="pi") %>% 
  mutate(chrom = "X")


Chromosome_X_fold_pi<-Chromosome_X_mean_pi_by_chr$mean_value[1]/Chromosome_X_mean_pi_by_chr$mean_value[2]
Chromosome_X_fold_pi
## 1.239609


















#################################################
####### Theta - Autosomes vs. Chromosome X ######
sub_region_raw<-read.table("../../data/05.07.21_cb_subregion.bed")

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

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



### Autosomes ##

Autosomes_theta_all<-window_diversity %>% 
  filter(chrom %in% c("I","II","III","IV","V")) %>% 
  filter(stat_type =="theta") %>% 
  mutate(mid = (window_start + window_stop) / 2) 



Autosomes_theta_with_region <- Autosomes_theta_all %>%
  inner_join(sub_region_df, by = "chrom") %>%
  filter(mid >= region_start, mid <= region_end)

Autosomes_theta_with_region <- Autosomes_theta_with_region %>%
  mutate(domain = case_when(
    subregion == "center"                 ~ "center",
    subregion %in% c("left_arm","right_arm") ~ "arm",
    TRUE                                  ~ NA_character_
  )) %>%
  filter(!is.na(domain))

Autosomes_mean_theta_by_chr <- Autosomes_theta_with_region %>%
  group_by(domain) %>%
  summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>% 
  mutate(type ="theta") %>% 
  mutate(chrom = "Autosomes")


Autosomes_fold_theta<-Autosomes_mean_theta_by_chr$mean_value[1]/Autosomes_mean_theta_by_chr$mean_value[2]
Autosomes_fold_theta
### 1.476646








### Chromosome_X ##

Chromosome_X_theta_all<-window_diversity %>% 
  filter(chrom %in% c("X")) %>% 
  filter(stat_type =="theta") %>% 
  mutate(mid = (window_start + window_stop) / 2) 



Chromosome_X_theta_with_region <- Chromosome_X_theta_all %>%
  inner_join(sub_region_df, by = "chrom") %>%
  filter(mid >= region_start, mid <= region_end)

Chromosome_X_theta_with_region <- Chromosome_X_theta_with_region %>%
  mutate(domain = case_when(
    subregion == "center"                 ~ "center",
    subregion %in% c("left_arm","right_arm") ~ "arm",
    TRUE                                  ~ NA_character_
  )) %>%
  filter(!is.na(domain))

Chromosome_X_mean_theta_by_chr <- Chromosome_X_theta_with_region %>%
  group_by(domain) %>%
  summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>% 
  mutate(type ="theta") %>% 
  mutate(chrom = "X")

Chromosome_X_fold_theta<-Chromosome_X_mean_theta_by_chr$mean_value[1]/Chromosome_X_mean_theta_by_chr$mean_value[2]
Chromosome_X_fold_theta
### 1.14169








#### #### #### #### #### 
#### output data #### 
#### #### #### #### #### 
all_table<-rbind(Autosomes_mean_pi_by_chr,
                 Chromosome_X_mean_pi_by_chr,
                 Autosomes_mean_theta_by_chr,
                 Chromosome_X_mean_theta_by_chr)

write.csv(all_table,
          "Autosomes_and_Chromosome_X_arm_vs_center.csv",
          quote = FALSE,row.names = FALSE)


Autosomes_fold_pi
# [1] 1.588569
Chromosome_X_fold_pi
# [1] 1.239459
Autosomes_fold_theta
# [1] 1.468403
Chromosome_X_fold_theta
# [1] 1.1412





