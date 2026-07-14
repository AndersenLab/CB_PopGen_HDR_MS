library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

ten_nigonis <- c("ECA2852","ECA2857","EG5268","JU1418","JU2617","JU4356","NIC2150","VX153","YR106","ZF1220")
aln <- readr::read_tsv("../../processed_data/gene_diversity/CBCN_nucmer_db_20250603.tsv",col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","CHROM","HIFI","STRAIN")) %>%
  dplyr::filter(STRAIN %in% ten_nigonis | STRAIN %in% c("ECA2666","ECA176")) %>%
  dplyr::filter(CHROM!="MtDNA") 

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

pn <- ggplot(aln %>% dplyr::filter(CHROM=="I",STRAIN %in% ten_nigonis)) + 
  geom_rect(data = region_rects2 %>% dplyr::filter(CHROM=="I"), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_rect(aes(xmin= S1/1e6, xmax=E1/1e6,ymin=-1,ymax=1),fill="black") + 
  facet_wrap(~STRAIN,ncol=1,strip.position = "top") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size=7,hjust = 1),
        panel.spacing.y = unit(0.01, "lines")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = region_colors) +
  xlab("Physical position (Mb)") +
  labs(fill="Chromosomal\ndomain")

pb <- ggplot(aln %>% dplyr::filter(CHROM=="I",!(STRAIN %in% ten_nigonis))) + 
  geom_rect(data = region_rects2 %>% dplyr::filter(CHROM=="I"), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_rect(aes(xmin= S1/1e6, xmax=E1/1e6,ymin=-1,ymax=1),fill="black") + 
  facet_wrap(~STRAIN,ncol=1,strip.position = "top") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size=7,hjust = 1),
        panel.spacing.y = unit(0.01, "lines"),
        axis.title.x=element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = region_colors) +
  xlab("Physical position (Mb)") +
  labs(fill="Chromosomal\ndomain") 


pcomp <- cowplot::plot_grid(pb,pn,ncol=1,axis = "v",align = "l",rel_heights = c(1,4),labels = c("a","b"))

ggsave(pcomp,filename = "../../figures/SF19_nigoni_accesibility.png",width = 7,height = 7,device = "png",units = "in",bg ="white",dpi = 900)

cre_aln <- readr::read_tsv('../../processed_data/tree_topology/CRE_files/CN_CE_transformed.tsv',col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","CHROM","HIFI"))

cre_plot <- ggplot() +
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_rect(data=cre_aln %>% dplyr::filter(CHROM!="MtDNA"),
            aes(xmin= S1/1e6, xmax=E1/1e6,ymin=-1,ymax=1),fill="black") + 
  facet_wrap(~CHROM,ncol=1,strip.position = "top",scales = "free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none",
        plot.title = element_text(hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size=7,hjust = 1),
        panel.spacing.y = unit(0.01, "lines")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = region_colors) +
  xlab("Physical position (Mb)") +
  labs(fill="Chromosomal\ndomain") 

ggsave(cre_plot,filename = "../../figures/SF20_remanei_accesibility.png",width = 7,height = 5,device = "png",units = "in",bg ="white",dpi = 900)
