library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

pi <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB/ALLCB_pi.txt")
theta <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB/ALLCB_watterson_theta.txt")
fst <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB/ALLCB_fst.txt")
tajD <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB/ALLCB_tajima_d.txt")
dxy <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB/ALLCB_dxy.txt")
dxy_silent <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB_4fold/ALLCB_4fold_dxy.txt")
global_pi <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB_nopop/ALLCB_nopop_pi.txt") %>% dplyr::mutate(pop="All")
global_theta <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB_nopop/ALLCB_nopop_watterson_theta.txt") %>% dplyr::mutate(pop="All")
global_tajD <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/pixy_ALLCB_nopop/ALLCB_nopop_tajima_d.txt") %>% dplyr::mutate(pop="All")
target_pops <- c("All","AD","KD","TD1","Temperate","TH","Tropical")
pop_sizes <- readr::read_tsv("../../processed_data/pixy_diversity_and_divergence/populations.txt",col_names = c("strains","pops")) %>%
  dplyr::group_by(pops) %>%
  dplyr::summarise(count=n())


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


pop_labeller <- c(
  "All" = "All (713)",
  pop_sizes %>%
    dplyr::filter(pops %in% target_pops) %>%
    dplyr::mutate(label = paste0(pops," (",count,")")) %>%
    dplyr::select(pops, label) %>%
    tibble::deframe()
)

piplot <- ggplot(rbind(global_pi,pi)%>% 
         dplyr::filter(pop %in% target_pops) %>%
         dplyr::rename(CHROM=chromosome) %>%
         dplyr::filter(CHROM!="MtDNA") %>%
         dplyr::mutate(
           pop = factor(
             pop,
             levels = c("All", "AD", "KD", "TD1", "Temperate", "TH", "Tropical")))) +
  facet_grid(
    pop ~ CHROM,
    scales = "free_x",
    labeller = ggplot2::labeller(pop = pop_labeller)
  ) +
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(aes(x=(window_pos_2-5e3)/1e6,y=avg_pi),size=0.05) + 
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 10),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=8)) +
        scale_x_continuous(expand = c(0,0)) +
  xlab("Physical position (Mb)") +
  ylab("Nucleotide diversity (π)")

ggsave(piplot, filename="../../figures/SF24_pi_est_by_group.png", width = 7, height = 7,device = "png",units = "in",bg = "white")


target_pi_stats <- rbind(global_pi,pi)%>% 
  dplyr::filter(pop %in% target_pops)
pi_dt <- as.data.table(target_pi_stats)
rect_dt <- as.data.table(region_rects)

rect_dt[, `:=`(
  domain_start = xmin * 1e6,
  domain_end   = xmax * 1e6
)]

chrom_max <- pi_dt[, .(chrom_end = max(window_pos_2, na.rm = TRUE)), by = chromosome]

rect_dt <- merge(
  rect_dt,
  chrom_max,
  by.x = "CHROM",
  by.y = "chromosome",
  all.x = TRUE
)

rect_dt[is.infinite(domain_end), domain_end := chrom_end]

# Prepare domains
domains_dt <- rect_dt[, .(
  chromosome = CHROM,
  region,
  domain_start,
  domain_end
)]

# Prepare windows
windows_dt <- copy(pi_dt)
windows_dt[, `:=`(
  query_start = window_pos_1,
  query_end   = window_pos_2
)]

setkey(domains_dt, chromosome, domain_start, domain_end)

classified <- foverlaps(
  x = windows_dt,
  y = domains_dt,
  by.x = c("chromosome", "query_start", "query_end"),
  by.y = c("chromosome", "domain_start", "domain_end"),
  type = "within",
  nomatch = 0L
)

window_domains <- as.data.frame(classified %>% dplyr::select(pop,chromosome,region,window_pos_1,window_pos_2))


thetaplot <- ggplot(rbind(global_theta,theta)%>% 
         dplyr::filter(pop %in% target_pops) %>%
         dplyr::rename(CHROM=chromosome) %>%
         dplyr::filter(CHROM!="MtDNA") %>%
         dplyr::mutate(
           pop = factor(
             pop,
             levels = c("All", "AD", "KD", "TD1", "Temperate", "TH", "Tropical")))) +
  facet_grid(
    pop ~ CHROM,
    scales = "free_x",
    labeller = ggplot2::labeller(pop = pop_labeller)
  ) +
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(aes(x=(window_pos_2-5e3)/1e6,y=avg_watterson_theta),size=0.05) + 
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 10),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=8)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.05)) +
  xlab("Physical position (Mb)") +
  ylab("Watterson's θ")


ggsave(thetaplot, filename="../../figures/SF25_theta_est_by_group.png", width = 7, height = 7,device = "png",units = "in",bg = "white")

target_pi_classified <- target_pi_stats %>%
  dplyr::left_join(window_domains, by=c("pop","chromosome", "window_pos_1", "window_pos_2")) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::rename(domain=region)


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
    dplyr::rename(CHROM = chromosome,ID=pop) %>%
    dplyr::filter(ID == id) %>%
    dplyr::filter(CHROM == chrom) %>%
    dplyr::mutate(x=(window_pos_1+window_pos_2)/2,
                  x_mb = x / x_divisor) %>%
    dplyr::filter(!is.na(avg_pi)) %>%
    dplyr::mutate(color=ifelse(avg_pi>=stat_thresh,"red","black"))
  
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
      aes(x = x_mb, y = avg_pi,
          color=color),
      size = point_size
    ) +
    geom_line(
      data = df,
      aes(x = x_mb, y = avg_pi)
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
targetstat <- target_pi_classified
ylab_tag = "π"
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
ylim_max = 0.025
SF12 <- cowplot::plot_grid(IL + theme(axis.title.x = element_blank()) + ylim(0,ylim_max),
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

ggsave(SF12,filename = "../../figures/SF12_localpi_tropical.png",width = 7,height = 6.5,device = "png",units = "in",bg ="white",dpi = 900)

summary_pi <- target_pi_classified %>%
  dplyr::mutate(
    Group = dplyr::recode(pop, All = "Global"),
    domain = dplyr::recode(domain, Tip = "Arm"),
    chrom_class = dplyr::if_else(chromosome == "X", "X", "Autosome")
  ) %>%
  dplyr::group_by(Group, chrom_class, domain) %>%
  dplyr::summarise(avg_pi = mean(avg_pi, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(col = paste(chrom_class, tolower(domain))) %>%
  dplyr::select(-chrom_class, -domain) %>%
  tidyr::pivot_wider(names_from = col, values_from = avg_pi) %>%
  dplyr::mutate(
    `fold change autosome arms/centers` = `Autosome arm` / `Autosome center`,
    `fold change X arms/centers` = `X arm` / `X center`,
    Stat = "\u03C0"
  ) %>%
  dplyr::rename(
    `Autosome arms` = `Autosome arm`,
    `Autosome centers` = `Autosome center`,
    `X arms` = `X arm`
  ) %>%
  dplyr::select(
    Group,
    `Autosome arms`,
    `Autosome centers`,
    `X arms`,
    `X center`,
    `fold change autosome arms/centers`,
    `fold change X arms/centers`,
    Stat
  ) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c("Global", "Tropical", "Temperate", "TH", "TD1", "KD", "AD")
    )
  ) %>%
  dplyr::arrange(Group) %>%
  dplyr::mutate(Group = as.character(Group)) %>%
  dplyr::left_join(pop_sizes %>% dplyr::rename(`Number of strains`=count),by=c("Group"="pops")) %>%
  dplyr::mutate(`Number of strains`=ifelse(is.na(`Number of strains`),713,`Number of strains`)) %>% 
  dplyr::mutate(Study="This study")


target_theta_stats <- rbind(global_theta,theta)%>% 
  dplyr::filter(pop %in% target_pops)

target_theta_classified <- target_theta_stats %>%
  dplyr::left_join(window_domains, by=c("pop","chromosome", "window_pos_1", "window_pos_2")) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::rename(domain=region)

summary_theta <- target_theta_classified %>%
  dplyr::mutate(
    Group = dplyr::recode(pop, All = "Global"),
    domain = dplyr::recode(domain, Tip = "Arm"),
    chrom_class = dplyr::if_else(chromosome == "X", "X", "Autosome")
  ) %>%
  dplyr::group_by(Group, chrom_class, domain) %>%
  dplyr::summarise(avg_watterson_theta  = mean(avg_watterson_theta , na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(col = paste(chrom_class, tolower(domain))) %>%
  dplyr::select(-chrom_class, -domain) %>%
  tidyr::pivot_wider(names_from = col, values_from = avg_watterson_theta ) %>%
  dplyr::mutate(
    `fold change autosome arms/centers` = `Autosome arm` / `Autosome center`,
    `fold change X arms/centers` = `X arm` / `X center`,
    Stat = "\u03B8W"
  ) %>%
  dplyr::rename(
    `Autosome arms` = `Autosome arm`,
    `Autosome centers` = `Autosome center`,
    `X arms` = `X arm`
  ) %>%
  dplyr::select(
    Group,
    `Autosome arms`,
    `Autosome centers`,
    `X arms`,
    `X center`,
    `fold change autosome arms/centers`,
    `fold change X arms/centers`,
    Stat
  ) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c("Global", "Tropical", "Temperate", "TH", "TD1", "KD", "AD")
    )
  ) %>%
  dplyr::arrange(Group) %>%
  dplyr::mutate(Group = as.character(Group)) %>%
  dplyr::left_join(pop_sizes %>% dplyr::rename(`Number of strains`=count),by=c("Group"="pops")) %>%
  dplyr::mutate(`Number of strains`=ifelse(is.na(`Number of strains`),713,`Number of strains`)) %>% 
  dplyr::mutate(Study="This study")


all_D <- rbind(tajD,global_tajD)

target_D_stats <- all_D %>% 
  dplyr::filter(pop %in% target_pops)

target_D_classified <- target_D_stats %>%
  dplyr::left_join(window_domains, by=c("pop","chromosome", "window_pos_1", "window_pos_2")) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::rename(domain=region)

summary_D <- target_D_classified %>%
  dplyr::mutate(
    Group = dplyr::recode(pop, All = "Global"),
    domain = dplyr::recode(domain, Tip = "Arm"),
    chrom_class = dplyr::if_else(chromosome == "X", "X", "Autosome")
  ) %>%
  dplyr::group_by(Group, chrom_class, domain) %>%
  dplyr::summarise(tajima_d = mean(tajima_d, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(col = paste(chrom_class, tolower(domain))) %>%
  dplyr::select(-chrom_class, -domain) %>%
  tidyr::pivot_wider(names_from = col, values_from = tajima_d ) %>%
  dplyr::mutate(
    `fold change autosome arms/centers` = "-",
    `fold change X arms/centers` = "-",
    Stat = "D_taj"
  ) %>%
  dplyr::rename(
    `Autosome arms` = `Autosome arm`,
    `Autosome centers` = `Autosome center`,
    `X arms` = `X arm`
  ) %>%
  dplyr::select(
    Group,
    `Autosome arms`,
    `Autosome centers`,
    `X arms`,
    `X center`,
    `fold change autosome arms/centers`,
    `fold change X arms/centers`,
    Stat
  ) %>%
  dplyr::mutate(
    Group = factor(
      Group,
      levels = c("Global", "Tropical", "Temperate", "TH", "TD1", "KD", "AD")
    )
  ) %>%
  dplyr::arrange(Group) %>%
  dplyr::mutate(Group = as.character(Group)) %>%
  dplyr::left_join(pop_sizes %>% dplyr::rename(`Number of strains`=count),by=c("Group"="pops")) %>%
  dplyr::mutate(`Number of strains`=ifelse(is.na(`Number of strains`),713,`Number of strains`)) %>% 
  dplyr::mutate(Study="This study")

all_stat_table <- rbind(summary_pi,summary_theta,summary_D)

write.table(all_stat_table,"../../supplementary_data/ST5_pi_theta_d_fold.tsv",quote = F,row.names = F,sep = "\t")

summary_genome <- target_pi_stats %>%
  dplyr::select(pop, chromosome, window_pos_1, window_pos_2, avg_pi) %>%
  dplyr::left_join(
    target_theta_stats %>%
      dplyr::select(pop, chromosome, window_pos_1, window_pos_2, avg_watterson_theta),
    by = c("pop", "chromosome", "window_pos_1", "window_pos_2")
  ) %>%
  dplyr::left_join(
    target_D_stats %>%
      dplyr::select(pop, chromosome, window_pos_1, window_pos_2, tajima_d),
    by = c("pop", "chromosome", "window_pos_1", "window_pos_2")
  ) %>%
  dplyr::group_by(pop) %>%
  dplyr::summarise(
    `π` = mean(avg_pi, na.rm = TRUE),
    `θW` = mean(avg_watterson_theta, na.rm = TRUE),
    `Tajima's D` = mean(tajima_d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    pop_sizes,
    by = c("pop" = "pops")
  ) %>%
  dplyr::mutate(
    Species = "C. briggsae",
    `Geographic region` = dplyr::recode(pop, All = "Global"),
    `Number of strains` = dplyr::if_else(pop == "All", 713L, count),
    Study = "This study",
    `Geographic region` = factor(
      `Geographic region`,
      levels = c("Global", "Tropical", "Temperate", "TH", "TD1", "KD", "AD")
    )
  ) %>%
  dplyr::arrange(`Geographic region`) %>%
  dplyr::mutate(`Geographic region` = as.character(`Geographic region`)) %>%
  dplyr::select(
    Species,
    `Geographic region`,
    `π`,
    `θW`,
    `Tajima's D`,
    `Number of strains`,
    Study
  )

write.table(summary_genome,"../../supplementary_data/ST4_pi_theta_d_mean.tsv",quote = F,row.names = F,sep = "\t")


tajdplot_A <- ggplot(all_D%>% 
         dplyr::filter(pop %in% c("All","Tropical")) %>%
         dplyr::rename(CHROM=chromosome) %>%
         dplyr::filter(CHROM!="MtDNA")%>%
         dplyr::mutate(
           pop = factor(
             pop,
             levels = c("All", "Tropical")))) +
  facet_grid(pop~CHROM,scales = "free_x")+
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(aes(x=(window_pos_2-5e3)/1e6,y=tajima_d),size=0.05) + 
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 10),
        legend.position = "none",
        panel.grid = element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y=element_text(size=10)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  xlab("Physical position (Mb)") +
  ylab("Tajima's D")

pop_order <- c("AD", "KD", "TD1", "TH", "Temperate", "Tropical", "All")
all_D2 <- all_D %>%
  dplyr::filter(pop %in% pop_order) %>%
  dplyr::mutate(
    pop = factor(pop, levels = pop_order),
    gid = as.numeric(pop)
  ) %>%
  dplyr::filter(!is.na(tajima_d)) %>%
  dplyr::filter(pop %in% target_pops) %>%
  dplyr::filter(chromosome!="MtDNA")

lut_lin <- all_D2 %>%
  dplyr::distinct(gid, pop) %>%
  dplyr::arrange(gid)

other_D_hm <- ggplot() +
  geom_rect(
    data = all_D2,
    aes(
      xmin = window_pos_1/1e6,
      xmax = window_pos_2/1e6,
      ymin = gid - 0.49,
      ymax = gid + 0.49,
      fill = tajima_d
    )
  ) +
  scale_fill_gradientn(
    colours = c("blue", "lightblue", "white", "pink", "red", "darkred"),
    values = scales::rescale(c(-2.5, -0.5, 0, 0.5, 2.5, 4), from = c(-2.5, 4)),
    limits = c(-2.5, 4),
    oob = scales::squish,
    breaks = c(-2, 0, 2, 4),
    name = "Tajima's D"
  )+
  ggh4x::facet_grid2(
    rows = vars(chromosome),
    scales = "free",
    space = "free_y"
  ) +
  scale_y_continuous(
    breaks = lut_lin$gid,
    labels = lut_lin$pop,
    expand = c(0,0),
    sec.axis = dup_axis(name = "")
  ) +
  labs(y = "Tajima's D", x = "Physical position (Mb)", fill = "Tajima's D") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        #strip.text.y = element_blank(),
        panel.border = element_rect(fill=NA),
        #strip.text = element_blank(),    
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.position = "bottom",        
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification.bottom = "right",
        legend.title.position = "left",
        legend.title = element_text(vjust = 1,size=8),
        #legend.title = element_blank(),
        plot.margin = margin(5, 5, 0, 10),
        legend.margin = margin(1, 1, 2, 5),
        axis.title.x=element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        legend.ticks = element_line(color = "black"),
        axis.title.y.right = element_text(size=10))



tajd_comb <- cowplot::plot_grid(tajdplot_A,other_D_hm,nrow=2,labels = c("a","b"),align = "v",axis = "lr",rel_heights = c(0.35,1))
tajd_comb <- cowplot::ggdraw(tajd_comb) +
  cowplot::draw_label("Physical position (Mb)",
                      x = 0.54, y = 0, vjust = -2, size = 10) 
  
ggsave(tajd_comb, filename="../../figures/SF26_tajD_by_group_hm.png", width = 7, height = 7.5,device = "png",units = "in",bg = "white")

tajdplot_all <- ggplot(all_D%>% 
                       dplyr::filter(pop %in% target_pops) %>%
                       dplyr::rename(CHROM=chromosome) %>%
                       dplyr::filter(CHROM!="MtDNA")%>%
                       dplyr::mutate(
                         pop = factor(
                           pop,
                           levels = c("All", "AD", "KD", "TD1", "Temperate", "TH", "Tropical")))) +
  facet_grid(pop~CHROM,scales = "free_x")+
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(aes(x=(window_pos_2-5e3)/1e6,y=tajima_d),size=0.05) + 
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 10),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=10)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(-5,5)) +
  xlab("Physical position (Mb)") +
  ylab("Tajima's D")

#ggsave(tajdplot_all, filename="../../figures/SFx_tajD_by_group_pt.png", width = 7, height = 7,device = "png",units = "in",bg = "white")

hdrs <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_allStrain_5kbclust_20260506.tsv") 

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

p1 <- ggplot(hdrs_ordered) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45)) + 
  facet_wrap(~CHROM,scales = 'free_x',nrow=1,) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  #axis.ticks.x=element_blank(),
  #axis.text.x=element_blank())  +
  ylab(paste0(length(unique(hdrs$STRAIN))," isotype strains")) +
  xlab("Physical position (Mb)") +
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 5),expand = c(0, 0))

global_stats <- rbind(global_pi %>% dplyr::mutate(stat_type="pi") %>% dplyr::select(pop,chromosome,window_start=window_pos_1,window_end=window_pos_2,stat=avg_pi,stat_type),
                      global_theta %>% dplyr::mutate(stat_type="theta") %>% dplyr::select(pop,chromosome,window_start=window_pos_1,window_end=window_pos_2,stat=avg_watterson_theta,stat_type)) %>%
  dplyr::rename(CHROM=chromosome) %>%
  dplyr::filter(CHROM!="MtDNA")

p2 <- ggplot(global_stats %>% 
               dplyr::mutate(stat_type = factor(ifelse(stat_type == "pi", "π", "Watterson's θ"), 
                                                           levels = c("π", "Watterson's θ")))) +
  facet_grid(
    stat_type ~ CHROM,
    scales = "free_x",
    labeller = ggplot2::labeller(pop = pop_labeller)
  ) +
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(aes(x=(window_end-5e3)/1e6,y=stat),size=0.05) + 
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        axis.title.x=element_blank(),
        strip.background.y = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(5, 5, 5, 10),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=8),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("Physical position (Mb)") +
  ylab("")
suppressWarnings(f3 <- cowplot::plot_grid(p2,p1,nrow=2,align="v",axis="lr",rel_heights = c(1,1.5), labels=c("a","b"))) # suppressing warnings of symbols in axis text (pi/theta) being substituted

ggsave(f3,filename = "../../figures/Figure3_nucdiversity_20260702.png",width = 7,height = 7.5,device = "png",units = "in",bg ="white",dpi = 900)

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
          dplyr::filter(check == FALSE & dplyr::coalesce(dplyr::lag(check), FALSE) == FALSE)
        
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


dxy_plot <- dxy %>%
  filter(pop1 %in% target_pops,
         pop2 %in% target_pops,
         chromosome != "MtDNA",
         pop1 == "AD",
         pop2 == "Tropical") %>%
  rename(CHROM = chromosome) %>%
  mutate(comp = ifelse(pop1 == "Tropical",
                       paste0("Tropical x ", pop2),
                       paste0("Tropical x ", pop1)))

chrom_means <- dxy_plot %>%
  group_by(CHROM) %>%
  summarize(mean_dxy = mean(avg_dxy, na.rm = TRUE),median_dxy = median(avg_dxy, na.rm = TRUE))

ggplot(dxy_plot) +
  facet_wrap(~CHROM, scales = "free_x", ncol = 1) +
  geom_rect(
    data = collapsed_tropical,
    aes(xmin = minStart/1e6, xmax = maxEnd/1e6, ymin = -Inf, ymax = Inf, fill = "grey60"),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_hline(
    data = chrom_means,
    aes(yintercept = mean_dxy),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  geom_line(
    aes(x = (window_pos_2 - 5e3) / 1e6,
        y = avg_dxy),
    linewidth = 0.2,
    color = "red"
  ) +
  scale_fill_manual(values = region_colors) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    strip.background.y = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 10),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 10)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.02)) +
  xlab("Physical position (Mb)") +
  ylab(expression(D[xy]))


annotate_overlaps <- function(combined_df, region_rects) {
  
  # Convert to data.table and prepare
  combined_dt <- as.data.table(combined_df %>%
                                 #filter(stat_type == "d") %>%
                                 mutate(window_pos_1 = window_pos_1 / 1e6,window_pos_2  = window_pos_2  / 1e6))
  
  region_dt <- as.data.table(region_rects)
  
  # Rename for foverlaps
  setnames(combined_dt, c("window_pos_1", "window_pos_2"), c("start", "end"))
  setnames(region_dt, c("xmin", "xmax"), c("start", "end"))
  
  # Ensure keys
  setkey(region_dt, CHROM, start, end)
  setkey(combined_dt, chromosome, start, end)
  
  # Overlap join
  annotated_dt <- foverlaps(region_dt, combined_dt, by.x = c("CHROM", "start", "end"), by.y = c("chromosome", "start", "end"), nomatch = 0, type = "any")
  return(annotated_dt)
}

tropical_d <- tajD  %>%
  dplyr::filter(pop=="Tropical")

annotated_by_lin <- annotate_overlaps(tropical_d, region_rects2) %>%
  dplyr::filter(pop=="Tropical") %>%
  dplyr::mutate(region_side = ifelse(side != "Center", paste0(side, "_", region), region)) %>%
  dplyr::filter(region != "Tip") %>%
  dplyr::mutate(
    region_side = factor(
      region_side,
      levels = c("Left_Arm", "Center", "Right_Arm"),
      labels = c("Left Arm", "Center", "Right Arm")
    )
  ) %>%
  dplyr::mutate(stat_type = "Tajima's D") %>%
  dplyr::mutate(x=start+((end-start)/2))

setDT(annotated_by_lin)
setDT(collapsed_tropical)

x <- copy(annotated_by_lin)[, `:=`(qstart = start, qend = end)]
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
ov[, `:=`(ov_len = pmax(0, pmin(qend, tend) - pmax(qstart, tstart)),
          int_len = pmax(0, qend - qstart))]
            
ov[, prop := fifelse(int_len > 0, ov_len / int_len, 0)]


ov_filt_best <- as.data.frame(ov) %>%
  dplyr::select(CHROM, x, start, end, prop) %>%
  dplyr::group_by(CHROM, x, start, end) %>%
  dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

anno_df_d <- as.data.frame(annotated_by_lin) %>%
  dplyr::left_join(ov_filt_best , by=c("CHROM","x","start","end")) %>%
  dplyr::mutate(hdstatus=ifelse(is.na(prop) | prop < 0.5, "non-HDR","HDR")) %>%
  dplyr::mutate(hdstatus = factor(hdstatus,
                                  levels = c("non-HDR","HDR"),
                                  labels = c("non-HDR","HDR")))

dlim=NA
tajdplot_poly <- ggplot(anno_df_d) +
  geom_freqpoly(aes(x = tajima_d, colour = hdstatus), binwidth = 0.05, linewidth = 1) +
  facet_grid(stat_type ~ region) +
  theme_bw(base_family = "Helvetica") +
  labs(color = "") +
  ylab("Genome span (Mb)") + xlab("") +
  coord_cartesian(xlim = c(-2.5, dlim)) +
  scale_y_continuous(labels = function(y) y / 100) +
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
  geom_jitter(aes(y = hdstatus, x = tajima_d, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = tajima_d),
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
  geom_jitter(aes(y = hdstatus, x = tajima_d, colour = hdstatus), alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = tajima_d),
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

tajdsumm <- cowplot::plot_grid(tajdplot_poly,
                               tajdbox,
                               nrow = 2,
                               rel_heights = c(2, 1),
                               align = "v",
                               axis = "lr")

obj <- readRDS("../../processed_data/tree_topology/aa_divergence/div_subplots.rds")

f6 <- cowplot::plot_grid(tajdsumm,cowplot::plot_grid(NULL,obj,nrow=1,rel_widths=c(0.04,1)),nrow=2,align = "v",axis = "lr",labels = c("a",""))

ggsave(f6, filename="../../figures/Figure6_D_Ks_Trop_CN.png", width = 7, height = 7,device = "png",units = "in",bg = "white")



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
    ) %>%
    dplyr::mutate(stat_type = "Dxy") %>%
    dplyr::mutate(x=start+((end-start)/2))
  
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

tropical_kd_dxy <- dxy  %>%
  dplyr::filter((pop1=="KD" | pop1 =="Tropical") & (pop2=="Tropical" | pop2=="KD"))

anno_dxy_kd <- process_dxy_annotation(tropical_kd_dxy, region_rects2)

tropical_kd_dxy_silent <- dxy_silent  %>%
  dplyr::filter((pop1=="KD" | pop1 =="Tropical") & (pop2=="Tropical" | pop2=="KD"))

anno_dxy_kd_silent <- process_dxy_annotation(tropical_kd_dxy_silent, region_rects2) %>%
  dplyr::mutate(hdstatus=ifelse(hdstatus=="HDR","HDR (syn)","non-HDR (syn)"))

dxyplot <- ggplot(rbind(anno_dxy_kd,anno_dxy_kd_silent) %>%
                    dplyr::filter(!is.na(avg_dxy))) +
  geom_freqpoly(aes(x = avg_dxy, colour = hdstatus), binwidth = 0.001, linewidth = 1) +
  facet_grid(stat_type ~ region) +
  theme_bw(base_family = "Helvetica") +
  labs(color = "") +
  ylab("Genome span (kb)") + xlab("") +
  coord_cartesian(xlim = c(0, 0.075)) +
  scale_y_continuous(labels = function(y) y * 10) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red","non-HDR (syn)" = "grey","HDR (syn)" = "pink")) +
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


dxybox <- ggplot(rbind(anno_dxy_kd,anno_dxy_kd_silent)) +
  geom_jitter(aes(y = hdstatus, x = avg_dxy, colour = hdstatus), alpha = 0.7) +
  geom_boxplot(aes(y = hdstatus, x = avg_dxy),
               colour = "grey40", 
               outliers = FALSE,
               width = 0.2, 
               outlier.shape = 16, 
               alpha = 0.6) +
  facet_grid(stat_type~region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab(expression(D[xy]))+
  coord_cartesian(xlim = c(0, 0.25)) +
  scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "red","non-HDR (syn)" = "grey","HDR (syn)" = "pink")) +
  theme(strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,5,5,5),
        text = element_text(family="Helvetica"),
        axis.title = element_text(size=10),
        axis.text = element_text(size=9),
        legend.text = element_text(size=9))

box_means <- rbind(anno_dxy_kd, anno_dxy_kd_silent) %>%
  group_by(stat_type, region, hdstatus) %>%
  summarise(
    mean_dxy = mean(avg_dxy, na.rm = TRUE),
    .groups = "drop"
  )

# dxybox_side <- ggplot(anno_dxy_kd) +
#   geom_jitter(aes(y = hdstatus, x = avg_dxy, colour = hdstatus), alpha = 0.2) +
#   geom_boxplot(aes(y = hdstatus, x = avg_dxy),
#                colour = "black", 
#                outliers = FALSE,
#                width = 0.2, 
#                outlier.shape = 16, 
#                alpha = 0.6) +
#   facet_wrap(~region,nrow=2,ncol=1,strip.position = "right") +
#   theme_bw(base_family = "Helvetica") +
#   ylab("") +
#   xlab(expression(D[xy]))+
#   #coord_cartesian(xlim = c(-2.5, dlim)) +
#   theme(legend.position = "none",
#         text = element_text(family="Helvetica"),
#         axis.title = element_text(size=10),
#         axis.text = element_text(size=9),
#         strip.text = element_text(size=10),
#         legend.text = element_text(size=9))

dxysumm <- cowplot::plot_grid(tajdplot_poly,
                              tajdbox,
                              dxyplot,
                              dxybox,
                              nrow = 4,
                              rel_heights = c(1.5, 1),
                              align = "v",
                              axis = "lr",
                              labels=c("a","","b","")) 




anno_dxy <- process_dxy_annotation(dxy, region_rects2)

anno_dxy_silent <- process_dxy_annotation(dxy_silent, region_rects2) %>%
  dplyr::mutate(hdstatus=ifelse(hdstatus=="HDR","HDR (syn)","non-HDR (syn)"))

all_dxy_wss <- rbind(anno_dxy,anno_dxy_silent) %>%
  dplyr::mutate(focalpop=ifelse(pop1=="Tropical",pop2,pop1)) %>%
  dplyr::filter(pop1=="Tropical" | pop2 =="Tropical") %>%
  dplyr::filter(!is.na(avg_dxy) & focalpop %in% target_pops)

my_comparisons <- list(
  c("non-HDR", "HDR"),             # nonsyn vs nonsyn
  c("non-HDR (syn)", "HDR (syn)")  # syn vs syn
)

dxybox_all <- ggplot(all_dxy_wss, aes(x = hdstatus, y = avg_dxy, colour = hdstatus)) +
  geom_jitter(alpha = 0.7, width = 0.3, height = 0,size=0.2) +
  geom_boxplot(
    colour = "grey40",
    outliers = FALSE,
    width = 0.3,
    alpha = 0.7
  ) +
  coord_flip(ylim = c(0, 0.1),)+
  facet_grid(focalpop ~ region) +
  theme_bw(base_family = "Helvetica") +
  ylab(expression(D[xy])) +
  xlab("") +
  scale_colour_manual(
    values = c(
      "non-HDR" = "black",
      "HDR" = "red",
      "non-HDR (syn)" = "grey",
      "HDR (syn)" = "pink"
    )
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5),
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    #axis.text.x = element_text(angle=60,hjust = 1),
    legend.text = element_text(size = 9)
  ) +
  scale_y_continuous(
    labels = function(x) sub("\\.?0+$", "", format(x, nsmall = 3))
  )

dxy_means2 <- all_dxy_wss %>%
  group_by(focalpop, region, hdstatus) %>%
  summarise(
    mean_dxy = mean(avg_dxy, na.rm = TRUE),
    n = dplyr::n(),
    sd = sd(avg_dxy, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

ST14 <- dxy_means2 %>%
  dplyr::mutate(`Group 1`="Tropical") %>%
  dplyr::rename(Domain=region) %>%
  tidyr::separate(hdstatus,into=c("region","sites"),sep = " \\(") %>%
  dplyr::mutate(sites=ifelse(is.na(sites),"All sites","Synonymous")) %>%
  dplyr::select(`Group 1`,`Group 2`=focalpop,Domain,Region=region,Sites=sites, `Mean Dxy`=mean_dxy) %>%
  tidyr::pivot_wider(
    names_from = Region,
    values_from = `Mean Dxy`) %>%
  dplyr::mutate(`Fold change (HDR / non-HDR)` = HDR / `non-HDR`) %>%
  dplyr::arrange(`Group 1`,`Group 2`,Sites,Domain) %>%
  dplyr::rename(`Mean Dxy (HDR)`=HDR,`Mean Dxy (non-HDR)`=`non-HDR`)

write.table(ST14,"../../supplementary_data/ST13_dxy_bypair.tsv",quote = F,row.names = F,sep = "\t")

ggsave(dxybox_all, filename="../../figures/SF27_dxy_to_Trop_all.png", width = 7, height = 7.7,device = "png",units = "in",bg = "white")
