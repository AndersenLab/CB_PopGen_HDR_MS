library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)


geo_colors <- c("Hawaii"="#66C2A5", 
                "Australia"="#FC8D62", 
                "Central America"="#8DA0CB",
                "South America"="#E78AC3", 
                "Africa"="#A6D854", 
                "Caribbean"="#FFD92F",
                "Taiwan" = "#E5C494", 
                "North America" = "#A65628", 
                "Europe" = "#E41A1C",
                "Asia" = "#399EB8", 
                "Pacific" = "#611EA1",
                "New Zealand" = "green",
                "Atlantic" = "purple", 
                "Oceania" ="#DB7779", 
                "Micronesia" = "#E7298A",
                "Indonesia" = "#7570B3", 
                "Malay Archipelago" = "#4110B3", 
                "C. America"="#8DA0CB",
                "S. America"="#E78AC3", 
                "N. America" = "#A65628", 
                "Globalized" = "gray30",
                "unknown" = 'grey')


all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/All_isotypes/chromosome_windows_diversity.csv",col_select = -1)

Temperate <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Temperate/chromosome_windows_diversity.csv",col_select = -1)
Tropical <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Tropical/chromosome_windows_diversity.csv",col_select = -1)
TH <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/TH/chromosome_windows_diversity.csv",col_select = -1)
TD1 <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/TD1/chromosome_windows_diversity.csv",col_select = -1)
AD <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/AD/chromosome_windows_diversity.csv",col_select = -1)
KD <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/KD/chromosome_windows_diversity.csv",col_select = -1)

Australia_all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Australia/chromosome_windows_diversity.csv",col_select = -1)
Caribbean_all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Caribbean/chromosome_windows_diversity.csv",col_select = -1)
Central_all <-readr::read_csv( "../../processed_data/diversity_and_divergence/pi_theta_d_geo/Central_America/chromosome_windows_diversity.csv",col_select = -1)
Taiwan_all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Taiwan/chromosome_windows_diversity.csv",col_select = -1)
Hawaii_all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Hawaii/chromosome_windows_diversity.csv",col_select = -1)
Pacific_all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Pacific/chromosome_windows_diversity.csv",col_select = -1)
SouthA_all <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/South_America/chromosome_windows_diversity.csv",col_select = -1)
Asia_all <-  readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_geo/Asia/chromosome_windows_diversity.csv",col_select = -1)

Taiwan_Trop <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Taiwan_Tropical/chromosome_windows_diversity.csv",col_select = -1)
Taiwan_TH <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Taiwan_TH/chromosome_windows_diversity.csv",col_select = -1)
Australia_Trop <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Australia_Tropical/chromosome_windows_diversity.csv",col_select = -1)
Asia_Trop <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Asia_Tropical/chromosome_windows_diversity.csv",col_select = -1)
Pacific_Trop <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Pacific_Tropical/chromosome_windows_diversity.csv",col_select = -1)
Hawaii_Trop <- readr::read_csv("../../processed_data/diversity_and_divergence/pi_theta_d_rg/Hawaii_Tropical/chromosome_windows_diversity.csv",col_select = -1)



Tropical$ID <- "Tropical"
Temperate$ID <- "Temperate"
AD$ID <- "AD"
KD$ID <- "KD"
TD1$ID <- "TD1"
TH$ID <- "TH"

Australia_all$ID <- "Australia"
Taiwan_all$ID <- "Taiwan"
Caribbean_all$ID <- "Caribbean"
Central_all$ID <- "Central America"
Asia_all$ID <- "Asia"
Hawaii_all$ID <- "Hawaii"
Pacific_all$ID <- "Pacific"
SouthA_all$ID <- "South America"


Taiwan_TH$ID <- "Taiwan (TH)"
Taiwan_Trop$ID <- "Taiwan (Tropical)"
Australia_Trop$ID <- "Australia (Tropical)"
Asia_Trop$ID <- "Asia (Tropical)"
Pacific_Trop$ID <- "Pacific (Tropical)"
Hawaii_Trop$ID <- "Hawaii (Tropical)"
all$ID <- "All isotypes"

# Combine all into one dataframe
combined_df_lin <- dplyr::bind_rows(
  all,
  Tropical,
  Temperate,
  KD,
  AD,
  TD1,
  TH
)

AD$ID <- "Australia (AD)"
KD$ID <- "Asia (KD)"
TD1$ID <- "Taiwan (TD1)"

combined_df_geo <- dplyr::bind_rows(
  Australia_all,
  Taiwan_all,
  Caribbean_all,
  Central_all,
  Asia_all,
  Hawaii_all,
  Pacific_all,
  SouthA_all)


Taiwan_all$ID <- "Taiwan (All)"
Australia_all$ID <- "Australia (All)"
Pacific_all$ID <- "Pacific (All)"
Asia_all$ID <- "Asia (All)"
Hawaii_all$ID <- "Hawaii (All)"

combined_df_geo2 <- dplyr::bind_rows(
  all,
  Australia_all,
  AD,
  Australia_Trop,
  Taiwan_all,
  TD1,
  Taiwan_TH,
  Taiwan_Trop,
  Caribbean_all,
  Central_all,
  Asia_all,
  Asia_Trop,
  KD,
  Hawaii_all,
  Hawaii_Trop,
  Pacific_all,
  Pacific_Trop,
  SouthA_all) %>%
  dplyr::mutate(group=ifelse(grepl("Hawaii",ID),2,
                             ifelse(grepl("Pacific",ID),3,
                                    ifelse(grepl("Taiwan",ID),4,
                                           ifelse(grepl("Asia",ID),5,
                                                  ifelse(grepl("Australia",ID),6,
                                                         ifelse(grepl("All",ID),0,1)))))))


#Read CSV and standardize column names/units
domains_raw <- readr::read_csv("../../processed_data/diversity_and_divergence/chromosome_domain_Cbriggsae.csv") %>%
  dplyr::rename(CHROM = chrom, start = left, end = right) %>%          
  dplyr::mutate(start = start * 1e3, end = end * 1e3)                 

#Create wide format with left and right arm domain coordinates
domains_wide <- domains_raw %>%
  dplyr::group_by(CHROM) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(arm = c("left", "right")) %>%                           
  tidyr::pivot_wider(                                                 
    names_from = arm,
    values_from = c(start, end),
    names_glue = "{arm}_{.value}") %>%
  dplyr::ungroup()

#Build rectangles for plotting chromosome regions: Tips, Arms, and Center
#Tip regions (from 0 to left_start and from right_end to Inf)
region_rects <- domains_wide %>%
  tidyr::pivot_longer(
    cols = c(left_start, right_end),
    names_to = "region_side",
    values_to = "x") %>%
  dplyr::mutate(
    region = "Tip",
    xmin = ifelse(region_side == "left_start", 0, x), 
    xmax = ifelse(region_side == "left_start", x, Inf)) %>%
  dplyr::select(CHROM, region, xmin, xmax) %>%
  #Add Arm regions
  dplyr::bind_rows(
    domains_wide %>%
      dplyr::mutate(region = "Arm") %>%
      dplyr::transmute(CHROM, region, xmin = left_start, xmax = left_end),
    domains_wide %>%
      dplyr::mutate(region = "Arm") %>%
      dplyr::transmute(CHROM, region, xmin = right_start, xmax = right_end),
    #Add Center region 
    domains_wide %>%
      dplyr::mutate(region = "Center") %>%
      dplyr::transmute(CHROM, region, xmin = left_end, xmax = right_start)) %>%
  #final formatting
  dplyr::mutate(
    ymin = -Inf,              
    ymax = Inf,
    xmin = xmin / 1e6,       
    xmax = xmax / 1e6)

region_rects2 <- region_rects %>% dplyr::arrange(region,CHROM,xmin) %>% dplyr::mutate(side=ifelse(row_number() %% 2 == 0,"Right","Left")) %>% dplyr::mutate(side=ifelse(region=="Center","Center",side))
region_colors <- c("Tip" = "#5E3C99", "Center" = "#FDB863", "Arm" = "#4393C3")

custom_order <- c("All isotypes",
                  "Caribbean",
                  "Central America",
                  "South America",
                  "Hawaii (All)",
                  "Hawaii (Tropical)",
                  "Pacific (All)",
                  "Pacific (Tropical)",
                  "Taiwan (All)",
                  "Taiwan (TD1)",
                  "Taiwan (TH)",
                  "Taiwan (Tropical)",
                  "Asia (All)",
                  "Asia (KD)",
                  "Asia (Tropical)",
                  "Australia (All)",
                  "Australia (AD)",
                  "Australia (Tropical)")  # whatever order you want

df_plot <- combined_df_geo2 %>%
  dplyr::filter(stat_type == "d") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(ID = factor(ID, levels = custom_order)) %>%
  dplyr::mutate(gid = as.integer(ID) * -1) %>%
  dplyr::ungroup()

df_plot_lin <- combined_df_lin %>%
  dplyr::filter(stat_type == "d") %>%
  dplyr::group_by(ID) %>%
  #mutate(ID = factor(ID, levels = custom_order)) %>%
  dplyr::mutate(gid = cur_group_id()) %>%
  dplyr::ungroup()

lut <- df_plot %>%
  dplyr::distinct(gid, ID) %>%
  dplyr::arrange(gid)

TDplot <- ggplot() +
  geom_rect(
    data = df_plot,
    aes(
      xmin = window_start/1e6,
      xmax = window_stop/1e6,
      ymin = gid - 0.49,
      ymax = gid + 0.49,
      fill = stat
    )
  ) +
  scale_fill_gradientn(
    colours = c("blue","lightblue", "white","pink", "red", "darkred"),  # blue → white → red → bright red
    values = scales::rescale(c(min(df_plot_lin$stat, na.rm = TRUE),
                               min(df_plot_lin$stat, na.rm = TRUE)/2,
                               0,
                               (min(df_plot_lin$stat, na.rm = TRUE)/2)*-1,
                               min(df_plot_lin$stat, na.rm = TRUE)*-1,
                               max(df_plot_lin$stat, na.rm = TRUE))),
    limits = c(min(df_plot_lin$stat, na.rm = TRUE), max(df_plot_lin$stat, na.rm = TRUE)),
    oob = scales::squish,  # keeps extreme values at ends
    name = ""
  ) +
  ggh4x::facet_grid2(
    rows = vars(group),
    cols = vars(chrom),
    scales = "free",
    space = "free_y" 
  ) +

  scale_y_continuous(
    breaks = lut$gid,
    labels = lut$ID,
    expand = c(0,0),
    sec.axis = dup_axis(name = "Tajima's D")
  ) +
  labs(y = "Tajima's D", x = "Physical position (Mb)", fill = "Tajima's D") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_family = "CyrHelvetica") +
  theme(axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.text = element_blank(),    
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.position = "bottom",        
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification.bottom = "right",
        legend.title.position = "top",
        legend.title = element_blank(),
        plot.margin = margin(0, 0, 0, 10),
        legend.margin = margin(0, 1, 0, 0),
        axis.title.x=element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.y.right = element_text(size=10),
        text=element_text(family="CyrHelvetica")) 

custom_order_lin <- c("All isotypes",
                  "Tropical",
                  "Temperate",
                  "TH",
                  "TD1",
                  "KD",
                  "AD")  # whatever order you want

df_plot_lin <- combined_df_lin %>%
  dplyr::filter(stat_type == "d") %>%
  dplyr::mutate(ID = factor(ID, levels = custom_order_lin)) %>%
  dplyr::mutate(gid = as.integer(ID) * -1) %>%
  dplyr::ungroup()

df_plot_lin_theta <- combined_df_lin %>%
  dplyr::filter(stat_type == "theta") %>%
  dplyr::mutate(ID = factor(ID, levels = custom_order_lin)) %>%
  dplyr::mutate(gid = as.integer(ID) * -1) %>%
  dplyr::ungroup()

df_plot_lin_pi <- combined_df_lin %>%
  dplyr::filter(stat_type == "pi") %>%
  dplyr::mutate(ID = factor(ID, levels = custom_order_lin)) %>%
  dplyr::mutate(gid = as.integer(ID) * -1) %>%
  dplyr::ungroup()


lut_lin <- df_plot_lin %>%
  dplyr::distinct(gid, ID) %>%
  dplyr::arrange(gid)

alliso_D_hm <- ggplot() +
  geom_rect(
    data = df_plot_lin %>% dplyr::filter(ID=="All isotypes"),
    aes(
      xmin = window_start/1e6,
      xmax = window_stop/1e6,
      ymin = gid - 0.49,
      ymax = gid + 0.49,
      fill = stat
    )
  ) +
  scale_fill_gradientn(
    colours = c("blue","lightblue", "white","pink", "red", "darkred"),  # blue → white → red → bright red
    values = scales::rescale(c(min(df_plot_lin$stat, na.rm = TRUE),
                               min(df_plot_lin$stat, na.rm = TRUE)/2,
                               0,
                               (min(df_plot_lin$stat, na.rm = TRUE)/2)*-1,
                               min(df_plot_lin$stat, na.rm = TRUE)*-1,
                               max(df_plot_lin$stat, na.rm = TRUE))),
    limits = c(min(df_plot_lin$stat, na.rm = TRUE), max(df_plot_lin$stat, na.rm = TRUE)),
    oob = scales::squish,  # keeps extreme values at ends
    name = ""
  ) +
  ggh4x::facet_grid2(
    cols = vars(chrom),
    scales = "free",
    space = "free_y" 
  ) +
  
  scale_y_continuous(
    breaks = lut_lin$gid,
    labels = lut_lin$ID,
    expand = c(0,0),
    sec.axis = dup_axis(name = "")
  ) +
  labs(y = "", x = "Physical position (Mb)", fill = "Tajima's D") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.text = element_blank(),    
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.position = "none",        
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification.bottom = "right",
        legend.title.position = "top",
        legend.title = element_blank(),
        plot.margin = margin(0, 0, 0, 10),
        legend.margin = margin(0, 1, 0, 0),
        axis.title.x=element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.ticks = element_line(color = "black"),
        axis.title.y.right = element_text(size=10)) 


trop_D_hm <- ggplot() +
  geom_rect(
    data = df_plot_lin %>% dplyr::filter(ID=="Tropical"),
    aes(
      xmin = window_start/1e6,
      xmax = window_stop/1e6,
      ymin = gid - 0.49,
      ymax = gid + 0.49,
      fill = stat
    )
  ) +
  scale_fill_gradientn(
    colours = c("blue","lightblue", "white","pink", "red", "darkred"),  # blue → white → red → bright red
    values = scales::rescale(c(min(df_plot_lin$stat, na.rm = TRUE),
                               min(df_plot_lin$stat, na.rm = TRUE)/2,
                               0,
                               (min(df_plot_lin$stat, na.rm = TRUE)/2)*-1,
                               min(df_plot_lin$stat, na.rm = TRUE)*-1,
                               max(df_plot_lin$stat, na.rm = TRUE))),
    limits = c(min(df_plot_lin$stat, na.rm = TRUE), max(df_plot_lin$stat, na.rm = TRUE)),
    oob = scales::squish,  # keeps extreme values at ends
    name = "Tajima's D"
  ) +
  ggh4x::facet_grid2(
    cols = vars(chrom),
    scales = "free",
    space = "free_y" 
  ) +
  
  scale_y_continuous(
    breaks = lut_lin$gid,
    labels = lut_lin$ID,
    expand = c(0,0),
    sec.axis = dup_axis(name = "")
  ) +
  labs(y = "", x = "Physical position (Mb)", fill = "Tajima's D") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.text = element_blank(),    
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.position = "none",        
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification.bottom = "right",
        legend.title.position = "top",
        legend.title = element_blank(),
        plot.margin = margin(0, 0, 0, 10),
        legend.margin = margin(0, 1, 0, 0),
        axis.title.x=element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.ticks = element_line(color = "black"),
        axis.title.y.right = element_text(size=10)) 


other_D_hm <- ggplot() +
  geom_rect(
    data = df_plot_lin,
    aes(
      xmin = window_start/1e6,
      xmax = window_stop/1e6,
      ymin = gid - 0.49,
      ymax = gid + 0.49,
      fill = stat
    )
  ) +
  scale_fill_gradientn(
    colours = c("blue","lightblue", "white","pink", "red", "darkred"),  # blue → white → red → bright red
    values = scales::rescale(c(min(df_plot_lin$stat, na.rm = TRUE),
                               min(df_plot_lin$stat, na.rm = TRUE)/2,
                               0,
                               (min(df_plot_lin$stat, na.rm = TRUE)/2)*-1,
                               min(df_plot_lin$stat, na.rm = TRUE)*-1,
                               max(df_plot_lin$stat, na.rm = TRUE))),
    limits = c(min(df_plot_lin$stat, na.rm = TRUE), max(df_plot_lin$stat, na.rm = TRUE)),
    oob = scales::squish,  # keeps extreme values at ends
    name = "Tajima's D  "
  ) +
  ggh4x::facet_grid2(
    cols = vars(chrom),
    scales = "free",
    space = "free_y" 
  ) +
  scale_y_continuous(
    breaks = lut_lin$gid,
    labels = lut_lin$ID,
    expand = c(0,0),
    sec.axis = dup_axis(name = "")
  ) +
  labs(y = "Tajima's D", x = "Physical position (Mb)", fill = "Tajima's D") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_rect(fill=NA),
        strip.text = element_blank(),    
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.position = "bottom",        
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification.bottom = "right",
        legend.title.position = "left",
        legend.title = element_text(vjust = 1,size=8),
        #legend.title = element_blank(),
        plot.margin = margin(15, 5, 5, 10),
        legend.margin = margin(2, 1, 2, 5),
        axis.title.x=element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.ticks = element_line(color = "black"),
        axis.title.y.right = element_text(size=10)) 


all_D <- ggplot() + 
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(data=df_plot_lin %>% 
               dplyr::filter(ID=="All isotypes") %>%
               dplyr::rename(CHROM=chrom),
             aes(x=x/1e6,y=stat),
             size=0.1) + 
  geom_smooth(data=df_plot_lin %>% 
                dplyr::filter(ID=="All isotypes") %>%
                dplyr::rename(CHROM=chrom),
              aes(x=x/1e6,y=stat),
              method="loess", se=FALSE, color="grey85", linewidth=0.6,span=0.3) + 
  facet_grid(ID~CHROM,scales = "free_x") +
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 10),
        legend.position = "none",
        #axis.title.y=element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=10)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill="Domain") +
  ylab("Tajima's D") +
  coord_cartesian(ylim=c(-2.5,4))

tropical_D <- ggplot() + 
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(data=df_plot_lin %>% 
               dplyr::filter(ID=="Tropical") %>%
               dplyr::rename(CHROM=chrom),
             aes(x=x/1e6,y=stat),
             size=0.1) + 
  geom_smooth(data=df_plot_lin %>% 
                dplyr::filter(ID=="Tropical") %>%
                dplyr::rename(CHROM=chrom),
              aes(x=x/1e6,y=stat),
              method="loess", se=FALSE, color="grey85", linewidth=0.6,span=0.3) + 
  facet_grid(ID~CHROM,scales = "free_x") +
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 1, 10),
        legend.position = "none",
        #axis.title.y=element_blank(),
        axis.text.y = element_text(size=9),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=10))+
        #strip.text=element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill="Domain") +
  ylab("Tajima's D") +
  coord_cartesian(ylim=c(-2.5,4))

EDF8 <- cowplot::plot_grid(all_D,tropical_D,other_D_hm,ncol=1,rel_heights = c(1,0.85,1.5),align = "v",axis="lr",labels=c("a","","b"))
EDF8 <- cowplot::ggdraw(EDF8) +
  cowplot::draw_label("Physical position (Mb)",
                      x = 0.54, y = 0, vjust = -2, size = 10) 

ggsave(EDF8,filename = "../../figures/EDF8_Dbylin_noAll.png",width = 7,height = 5.5,device = "png",units = "in",bg ="white",dpi = 900)


plot_lin_pi_region <- function(region_rects,
                               stats_df,
                               chrom = "III",
                               id = "Tropical",
                               side = "R",
                               x_lim_mb = c(10, 14),
                               region_colors = NULL,
                               point_size = 0.5,
                               rect_alpha = 0.5,
                               x_divisor = 1e6,
                               show_legend = FALSE,
                               y_lab = "π",
                               title = NULL,
                               stat_thresh = 0.002) {
  # Prep data
  rr <- region_rects %>%
    dplyr::filter(CHROM == chrom)
  
  df <- stats_df %>%
    dplyr::filter(ID == id) %>%
    dplyr::rename(CHROM = chrom) %>%
    dplyr::filter(CHROM == chrom) %>%
    dplyr::mutate(x_mb = x / x_divisor) %>%
    dplyr::filter(!is.na(stat)) %>%
    dplyr::mutate(color=ifelse(stat>=stat_thresh,"red","black"))
  
  # Default title
  if (is.null(title)) {
    title <- sprintf("%s - %s", chrom, side)
  }
  
  p <- ggplot() +
    geom_hline(yintercept = stat_thresh,linetype="dashed",color="grey")+
    geom_rect(
      data = rr,
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax, fill = region),
      inherit.aes = FALSE, alpha = rect_alpha
    ) +
    geom_point(
      data = df,
      aes(x = x_mb, y = stat,
          color=color),
      size = point_size
    ) +
    geom_line(
      data = df,
      aes(x = x_mb, y = stat)
    ) +
    facet_wrap(~CHROM, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = region_colors) +
    theme_bw() +
    theme(
      panel.border = element_rect(fill = NA),
      strip.background.y = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      plot.margin = margin(1,10, 2, 2),
      plot.title = element_text(size=9,vjust = 1,margin = margin(b = 1)),
      legend.position = if (show_legend) "right" else "none",
      axis.text.y = element_text(size = 9),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(fill = "Domain") +
    coord_cartesian(xlim = x_lim_mb) +
    ggtitle(title) +
    xlab("Physical position (Mb)") +
    ylab(y_lab) +
    scale_color_identity()
  
  return(p)
}

isogroup = "Tropical"
targetstat <- df_plot_lin_pi
ylab_tag = "π"
#ylab_tag = "θ"
stat_thresh = 0.002

IR <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "I",
  id = isogroup,
  side = "R",
  x_lim_mb = c(9.2, 13.2),
  region_colors = region_colors
)

IL <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "I",
  id = isogroup,
  x_lim_mb = c(1.2, 5.2),
  side = "L",
  region_colors = region_colors
)


IIR <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "II",
  id = isogroup,
  side = "R",
  x_lim_mb = c(12, 16),
  region_colors = region_colors
)

IIL <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "II",
  id = isogroup,
  x_lim_mb = c(1.4, 5.4),
  side = "L",
  region_colors = region_colors
)

IIIR <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "III",
  id = isogroup,
  side = "R",
  x_lim_mb = c(9.8, 13.8),
  region_colors = region_colors
)

IIIL <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "III",
  id = isogroup,
  x_lim_mb = c(1.9, 5.9),
  side = "L",
  region_colors = region_colors
)


IVR <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "IV",
  id = isogroup,
  side = "R",
  x_lim_mb = c(11.4, 15.4),
  region_colors = region_colors
)

IVL <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "IV",
  id = isogroup,
  x_lim_mb = c(2.2, 6.2),
  side = "L",
  region_colors = region_colors
)

VR <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "V",
  id = isogroup,
  side = "R",
  x_lim_mb = c(13.5, 17.5),
  region_colors = region_colors
)

VL <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "V",
  id = isogroup,
  x_lim_mb = c(4, 8),
  side = "L",
  region_colors = region_colors
)

XR <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "X",
  id = isogroup,
  side = "R",
  x_lim_mb = c(16.5, 20.5),
  region_colors = region_colors
)

XL <- plot_lin_pi_region(
  region_rects = region_rects2,
  stats_df = targetstat,
  y_lab = ylab_tag,
  stat_thresh = stat_thresh,
  chrom = "X",
  id = isogroup,
  x_lim_mb = c(5, 9),
  side = "L",
  region_colors = region_colors
)
ylim_max = 0.013
SF8 <- cowplot::plot_grid(IL + theme(axis.title.x = element_blank()) + ylim(0,ylim_max),
                   IR + theme(axis.title = element_blank()) + ylim(0,ylim_max),
                   IIL + theme(axis.title.x = element_blank()) + ylim(0,ylim_max),
                   IIR + theme(axis.title = element_blank()) + ylim(0,ylim_max),
                   IIIL + theme(axis.title.x = element_blank()) + ylim(0,ylim_max),
                   IIIR + theme(axis.title = element_blank()) + ylim(0,ylim_max),
                   IVL + theme(axis.title.x = element_blank()) + ylim(0,ylim_max),
                   IVR + theme(axis.title = element_blank()) + ylim(0,ylim_max),
                   VL + theme(axis.title.x = element_blank())+ ylim(0,ylim_max),
                   VR + theme(axis.title = element_blank()) + ylim(0,ylim_max),
                   XL + ylim(0,ylim_max),
                   XR + theme(axis.title.y = element_blank()) + ylim(0,ylim_max),
                   align = "v",axis="lr",ncol=2,rel_heights = c(1,1,1,1,1,1.2))

ggsave(SF8,filename = "../../figures/SF8_localpi_tropical.png",width = 7,height = 6.5,device = "png",units = "in",bg ="white",dpi = 900)


pitheta_plot <- ggplot() + 
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(data=all %>% 
               dplyr::filter(stat_type!="d" & stat_type!="mis") %>%
               dplyr::rename(CHROM=chrom) %>%
               dplyr::mutate(stat_type = factor(ifelse(stat_type == "pi", "π", "Watterson's θ"), 
                                                levels = c("π", "Watterson's θ"))),
             aes(x=x/1e6,y=stat),
             size=0.05) +
  geom_smooth(data=all %>% 
                dplyr::filter(stat_type!="d" & stat_type!="mis") %>%
                dplyr::rename(CHROM=chrom) %>%
                dplyr::mutate(stat_type = factor(ifelse(stat_type == "pi", "π", "Watterson's θ"), 
                                                 levels = c("π", "Watterson's θ"))),
              aes(x=x/1e6,y=stat),
              method="loess", se=FALSE, color="grey85", linewidth=1,span=0.3) + 
  scale_fill_manual(values = region_colors) +
  facet_grid(stat_type~CHROM,scales = "free_x") +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(0, 0, 5, 10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=10)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill="Domain")



F3_legacy <- cowplot::plot_grid(pitheta_plot,TDplot,nrow=2,ncol=1,rel_heights = c(1,1.3),align = "v",axis = "lr")
F3_legacy <- cowplot::ggdraw(F3_legacy) +
  cowplot::draw_label("Physical position (Mb)",
                      x = 0.54, y = 0, vjust = -1, size = 10)

hdrs <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_allStrain_5kbclust_20250930.tsv") 

hdrs_ordered <- hdrs %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() 

# hdrs_unt_rg_ordered <- hdrs_unt_rg %>% 
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(ncalls=n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::arrange(source,desc(ncalls),STRAIN,CHROM,minStart) %>%
#   dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
#   dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(ystrain=cur_group_id()) %>%
#   dplyr::ungroup()
# 
# x_breaks_by_chrom <- bins_wFreq %>%
#   group_by(CHROM) %>%
#   summarise(x_breaks = list(seq(
#     floor(min((start) / 1e6)),
#     ceiling(max((end) / 1e6)),
#     by = 5
#   )))

p1 <- ggplot(hdrs_ordered) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45)) + 
  facet_wrap(~CHROM,scales = 'free_x',nrow=1) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank())  +
  ylab("695 Isotype strains") +
  xlab("Physical position (Mb)") +
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 5),expand = c(0, 0))


F3_alt <- cowplot::plot_grid(pitheta_plot,p1,ncol=1,align="v",axis="lr",rel_heights = c(1,1.5), labels=c("a","b"))

ggsave(F3_alt,filename = "../../figures/Figure3_ALT_nucdiversity_20251002.png",width = 7,height = 7.5,device = "png",units = "in",bg ="white",dpi = 900)


getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      #print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>% 
        dplyr::arrange(CHROM,minStart) %>%
        dplyr::mutate(check=ifelse(lead(minStart) <= maxEnd,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
      if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
        print("NO MORE INTERSECTS")
        k=0
      } else {
        
        temp <- checkIntersect %>%
          dplyr::mutate(gid=data.table::rleid(check)) %>%
          dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))
        
        collapse <- temp %>%
          dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart=min(minStart)) %>%
          dplyr::mutate(newEnd=max(maxEnd)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(minStart=newStart,maxEnd=newEnd) %>%
          dplyr::select(-newEnd,-newStart)
        
        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print("WRITING TO LIST")
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

collapsed_tropical <- plyr::ldply(getRegFreq(hdrs %>% 
                                               dplyr::filter(source=="QX1410") %>% 
                                               dplyr::group_split(CHROM)),data.frame) %>% 
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::select(-STRAIN)

cov_wg <- sum(collapsed_tropical$divSize) / 106184000
wg_tot <- nrow(collapsed_tropical)
wg_mean <- mean(collapsed_tropical$divSize) / 1e3
wg_min <- min(collapsed_tropical$divSize) / 1e3
wg_max <- max(collapsed_tropical$divSize) / 1e3

bins <- readr::read_tsv("../../processed_data/HDRs/QX1410_genomic_windows.1kb.bed",col_names = c("CHROM","binStart","binEnd")) 

bins_dt <- as.data.table(bins)
setnames(bins_dt, c("binStart", "binEnd"), c("start", "end"))
bins_dt[, id := .I]  # optional: keep track of bins

hdrs_dt <- as.data.table(hdrs %>% dplyr::filter(source=="QX1410"))
setnames(hdrs_dt, c("minStart", "maxEnd"), c("start", "end"))

setkey(bins_dt, CHROM, start, end)
setkey(hdrs_dt, CHROM, start, end)

overlaps <- foverlaps(hdrs_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(STRAIN)), by = .(CHROM, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("CHROM", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq=n_strains/length(unique((hdrs %>% dplyr::filter(source=="QX1410"))$STRAIN)))

x_breaks_by_chrom <- bins_wFreq %>%
  group_by(CHROM) %>%
  summarise(x_breaks = list(seq(
    floor(min((start) / 1e6)),
    ceiling(max((end) / 1e6)),
    by = 5
  )))

freq_plot <- ggplot(bins_wFreq) +
  geom_point(aes(x=(start+500)/1e6,y=freq),size=0.2,stroke = 0) +
  facet_wrap(~CHROM,scales = 'free_x',nrow=1) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        strip.text = element_blank(),
        plot.margin=margin(t=0))  +
  ylab("Tropical group\nHDR frequency")+
  xlab("Physical position (Mb)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x+(500/1e6))), by = 5),expand = c(0, 0))
freq_plot


hdrs_panel <- ggplot(hdrs_ordered) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45)) + 
  facet_wrap(~CHROM,scales = 'free_x',nrow=1) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank())  +
  ylab("695 Isotype strains") +
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 5),expand = c(0, 0))

F3_alt_wfreq <- cowplot::plot_grid(pitheta_plot,hdrs_panel,freq_plot,ncol=1,align="v",axis="lr",rel_heights = c(1,1.5,0.5), labels=c("a","b",""))

#ggsave(F3_alt_wfreq,filename = "../../figures/Figure3_ALT_nucdiversity_wfreq_20251002.png",width = 7,height = 7.5,device = "png",units = "in",bg ="white",dpi = 900)

annotate_overlaps <- function(combined_df, region_rects) {

  # Convert to data.table and prepare
  combined_dt <- as.data.table(combined_df %>%
                                 #filter(stat_type == "d") %>%
                                 mutate(window_start = window_start / 1e6,window_stop  = window_stop  / 1e6))

  region_dt <- as.data.table(region_rects)

  # Rename for foverlaps
  setnames(combined_dt, c("window_start", "window_stop"), c("start", "end"))
  setnames(region_dt, c("xmin", "xmax"), c("start", "end"))

  # Ensure keys
  setkey(region_dt, CHROM, start, end)
  setkey(combined_dt, chrom, start, end)

  # Overlap join
  annotated_dt <- foverlaps(region_dt, combined_dt, by.x = c("CHROM", "start", "end"), by.y = c("chrom", "start", "end"), nomatch = 0, type = "any")
  return(annotated_dt)
}

annotated_by_lin <- annotate_overlaps(df_plot_lin, region_rects2) %>%
  dplyr::filter(ID=="Tropical") %>%
  dplyr::mutate(region_side = ifelse(side != "Center", paste0(side, "_", region), region)) %>%
  dplyr::filter(region != "Tip") %>%
  dplyr::mutate(
    region_side = factor(
      region_side,
      levels = c("Left_Arm", "Center", "Right_Arm"),
      labels = c("Left Arm", "Center", "Right Arm")
    )
  )

df_d <- annotated_by_lin %>%
      dplyr::filter(stat_type == "d") %>%
      dplyr::mutate(stat_type = "Tajima's D")

setDT(df_d)
setDT(collapsed_tropical)

# make copies so we don't clobber original columns
x <- copy(df_d)[, `:=`(qstart = start, qend = end)]
y <- copy(collapsed_tropical)[
  , `:=`(tstart = minStart / 1e6, tend = maxEnd / 1e6)
][, .(CHROM, tstart, tend)]

# add a unique row id to df_d to aggregate after overlaps
x[, rowid := .I]

# set keys for foverlaps
setkey(x, CHROM, qstart, qend)
setkey(y, CHROM, tstart, tend)

# all overlaps on matching CHROM
ov <- foverlaps(
  x = x, y = y,
  by.x = c("CHROM", "qstart", "qend"),
  by.y = c("CHROM", "tstart", "tend"),
  type = "any", nomatch = 0L
)

# compute overlap proportion relative to the df_d interval
ov[, `:=`(
  ov_len = pmax(0, pmin(qend, tend) - pmax(qstart, tstart)),
  int_len = pmax(0, qend - qstart)
)]
ov[, prop := fifelse(int_len > 0, ov_len / int_len, 0)]


ov_filt_best <- as.data.frame(ov) %>%
  dplyr::select(CHROM, x, start, end, prop) %>%
  dplyr::group_by(CHROM, x, start, end) %>%
  dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

anno_df_d <- as.data.frame(df_d) %>%
  dplyr::left_join(ov_filt_best , by=c("CHROM","x","start","end")) %>%
  dplyr::mutate(hdstatus=ifelse(is.na(prop) | prop < 0.5, "non-HDR","HDR")) %>%
  dplyr::mutate(hdstatus = factor(hdstatus,
                                  levels = c("non-HDR","HDR"),
                                  labels = c("non-HDR","HDR")))




dlim=NA

tajdplot <- ggplot(anno_df_d) +
  geom_freqpoly(aes(x = stat, colour = hdstatus), binwidth = 0.05, linewidth = 1) +
  facet_grid(stat_type ~ region) +
  theme_bw(base_family = "Helvetica") +
  labs(color = "") +
  ylab("Genome span (kb)") + xlab("") +
  coord_cartesian(xlim = c(-2.5, dlim)) +
  scale_y_continuous(labels = function(y) y * 10) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text.y = element_blank(),
        legend.position = "inside",
        legend.justification.inside = c(0.99, 0.95),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5,5,0,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        strip.text = element_text(size=10),
        legend.text = element_text(size=9))

tajdbox <- ggplot(anno_df_d) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_wrap(stat_type~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab("Tajima's D") +
  coord_cartesian(xlim = c(-2.5, dlim)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

tajdbox_side <- ggplot(anno_df_d) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "black", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_wrap(~region,nrow=2,ncol=1,strip.position = "right") +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab("Tajima's D") +
  coord_cartesian(xlim = c(-2.5, dlim)) +
  theme(legend.position = "none")

tajdsumm <- cowplot::plot_grid(tajdplot,
                               tajdbox,
                               nrow = 2,
                               rel_heights = c(2, 1),
                               align = "v",
                               axis = "lr")  #   tajdsumm <- cowplot::plot_grid(tajdplot,

dxy_tro_ad <- readr::read_delim("../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_AD__dxy.txt") %>%
  dplyr::rename(chrom=chromosome,window_start=window_pos_1,window_stop=window_pos_2,stat=avg_dxy)  %>%
  dplyr::mutate(comp_lin="AD")

dxy_tro_kd <- readr::read_delim("../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_KD__dxy.txt") %>%
  dplyr::rename(chrom=chromosome,window_start=window_pos_1,window_stop=window_pos_2,stat=avg_dxy)  %>%
  dplyr::mutate(comp_lin="KD")

dxy_tro_td <- readr::read_delim("../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_TD1__dxy.txt") %>%
  dplyr::rename(chrom=chromosome,window_start=window_pos_1,window_stop=window_pos_2,stat=avg_dxy)  %>%
  dplyr::mutate(comp_lin="TD1")

dxy_tro_tmp <- readr::read_delim("../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_Temperate__dxy.txt") %>%
  dplyr::rename(chrom=chromosome,window_start=window_pos_1,window_stop=window_pos_2,stat=avg_dxy)  %>%
  dplyr::mutate(comp_lin="Temperate")

dxy_tro_th <- readr::read_delim("../../processed_data/diversity_and_divergence/Dxy_Tropical/Tropical_TH__dxy.txt") %>%
  dplyr::rename(chrom=chromosome,window_start=window_pos_1,window_stop=window_pos_2,stat=avg_dxy) %>%
  dplyr::mutate(comp_lin="TH")



process_dxy_annotation <- function(dxy_tro_lin, region_rects2) {
  # Annotate overlaps and prepare annotated data
  annotated_dxy_lin <- annotate_overlaps(dxy_tro_lin, region_rects2) %>%
    dplyr::mutate(region_side = ifelse(side != "Center", paste0(side, "_", region), region)) %>%
    dplyr::filter(region != "Tip") %>%
    dplyr::mutate(
      region_side = factor(
        region_side,
        levels = c("Left_Arm", "Center", "Right_Arm"),
        labels = c("Left Arm", "Center", "Right Arm")
      )
    )
  
  # Convert to data.table for overlap operations
  setDT(annotated_dxy_lin)
  x <- copy(annotated_dxy_lin)[, `:=`(qstart = start, qend = end)]
  setkey(x, CHROM, qstart, qend)
  
  # Perform overlap with y
  ov_dxy <- foverlaps(
    x = x, y = y, #y was previously defined as a dt of collapsed tropical df
    by.x = c("CHROM", "qstart", "qend"),
    by.y = c("CHROM", "tstart", "tend"),
    type = "any", nomatch = 0L
  )
  
  # Compute overlap proportion
  ov_dxy[, `:=`(
    ov_len = pmax(0, pmin(qend, tend) - pmax(qstart, tstart)),
    int_len = pmax(0, qend - qstart)
  )]
  ov_dxy[, prop := fifelse(int_len > 0, ov_len / int_len, 0)]
  
  # Select best overlaps
  ov_dxy_best <- as.data.frame(ov_dxy) %>%
    dplyr::select(CHROM, start, end, prop) %>%
    dplyr::group_by(CHROM, start, end) %>%
    dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Join with annotated data and assign HDR status
  anno_dxy <- as.data.frame(annotated_dxy_lin) %>%
    dplyr::left_join(ov_dxy_best, by = c("CHROM", "start", "end")) %>%
    dplyr::mutate(
      hdstatus = ifelse(is.na(prop) | prop < 0.5, "non-HDR", "HDR"),
      hdstatus = factor(hdstatus, levels = c("non-HDR", "HDR")),
      stat_type = "Dxy"
    )
  
  return(anno_dxy)
}

anno_dxy_ad <- process_dxy_annotation(dxy_tro_ad, region_rects2)

dxyplot <- ggplot(anno_dxy_ad) +
  geom_freqpoly(aes(x = stat, colour = hdstatus), binwidth = 0.01, linewidth = 1) +
  facet_grid(stat_type ~ region) +
  theme_bw(base_family = "Helvetica") +
  labs(color = "") +
  ylab("Genome span (kb)") + xlab("") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_y_continuous(labels = function(y) y * 10) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text.y = element_blank(),
        legend.position = "inside",
        legend.justification.inside = c(0.99, 0.95),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(5,5,0,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        strip.text = element_text(size=10),
        legend.text = element_text(size=9))


dxybox <- ggplot(anno_dxy_ad) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(stat_type~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab(expression(D[xy]))+
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

dxybox_side <- ggplot(anno_dxy_ad) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "black", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_wrap(~region,nrow=2,ncol=1,strip.position = "right") +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab(expression(D[xy]))+
  #coord_cartesian(xlim = c(-2.5, dlim)) +
  theme(legend.position = "none",
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        strip.text = element_text(size=10),
        legend.text = element_text(size=9))

dxysumm <- cowplot::plot_grid(tajdplot,
                              tajdbox,
                              dxyplot,
                              dxybox,
                               nrow = 4,
                               rel_heights = c(1.5, 1),
                               align = "v",
                               axis = "lr",
                              labels=c("a","","b","")) 
# 




ggsave(dxysumm,filename = "../../figures/Figure5_TajD_Dxy_20251007.png",width = 7,height = 7.5,device = "png",units = "in",bg ="white",dpi = 900)


anno_dxy_ad <- process_dxy_annotation(dxy_tro_ad, region_rects2)
dxybox_ad <- ggplot(anno_dxy_ad) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(comp_lin~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab("")+
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(5,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

anno_dxy_kd <- process_dxy_annotation(dxy_tro_kd, region_rects2)
dxybox_kd <- ggplot(anno_dxy_kd) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(comp_lin~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab("")+
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))


anno_dxy_th <- process_dxy_annotation(dxy_tro_th, region_rects2)
dxybox_th <- ggplot(anno_dxy_th) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(comp_lin~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab("")+
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

anno_dxy_td <- process_dxy_annotation(dxy_tro_td, region_rects2)
dxybox_td <- ggplot(anno_dxy_td) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(comp_lin~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab("")+
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

anno_dxy_tmp <- process_dxy_annotation(dxy_tro_tmp, region_rects2)
dxybox_tmp <- ggplot(anno_dxy_tmp) +
  geom_jitter(aes(y = hdstatus, x = stat, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(comp_lin~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab(expression(D[xy]))+
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red")) +
  theme(strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

dxy_boxes <- cowplot::plot_grid(dxybox_ad,dxybox_kd,dxybox_td,dxybox_th,dxybox_tmp,nrow=5,rel_heights = c(1,0.85,0.85,0.85,0.85))

ggsave(dxy_boxes,filename = "../../figures/EDF9_Dxy_allcomp_20251007.png",width = 7,height = 7.5,device = "png",units = "in",bg ="white",dpi = 900)
