library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

tropvtrop <- readr::read_tsv("../../processed_data/HDRs/Tropical/counts/ECA2645.variant_counts.tsv",col_names=c("chrom","start","end","count")) %>%
  dplyr::filter(chrom %in% c("I","III","V"))
tropvkd <- readr::read_tsv("../../processed_data/HDRs/all_strain/counts/ECA2666.variant_counts.tsv",col_names=c("chrom","start","end","count")) %>%
  dplyr::filter(chrom %in% c("I","III","V"))
kdvkd <- readr::read_tsv("../../processed_data/HDRs/Other_RG/counts/ECA2670/ECA2666.variant_counts.tsv",col_names=c("chrom","start","end","count")) %>%
  #dplyr::filter(!grepl("ptg|MtDNA",chrom)) %>%
  dplyr::filter(chrom %in% c("I","III","V"))

tropvtrop <- readr::read_tsv("../../processed_data/HDRs/Tropical/counts/ECA2645.variant_counts.tsv",col_names=c("chrom","start","end","count"))
tropvkd <- readr::read_tsv("../../processed_data/HDRs/all_strain/counts/ECA2666.variant_counts.tsv",col_names=c("chrom","start","end","count")) 
kdvkd <- readr::read_tsv("../../processed_data/HDRs/Other_RG/counts/ECA2670/ECA2666.variant_counts.tsv",col_names=c("chrom","start","end","count")) %>%
  dplyr::filter(!grepl("ptg|MtDNA",chrom))
  #dplyr::filter(chrom %in% c("I","III","V"))


p1 <- ggplot(tropvtrop) + 
  geom_point(aes(x=(start+500)/1e6,y=count),size=0.1) + 
  facet_wrap(~chrom,nrow=1,scales="free_x") +
  ylab("SNVs / kb") +
  xlab("") +
  #ggtitle("Tropical vs. Tropical") +
  theme_bw()  +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,100),expand = c(0.01,0))

p2 <- ggplot(tropvkd) + 
  geom_point(aes(x=(start+500)/1e6,y=count),size=0.1) + 
  facet_wrap(~chrom,nrow=1,scales="free_x") +
  ylab("SNVs / kb") +
  xlab("") +
  #ggtitle("Tropical vs. KD ")+
  theme_bw()  +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,100),expand = c(0.01,0))

p3 <- ggplot(kdvkd) + 
  geom_point(aes(x=(start+500)/1e6,y=count),size=0.1) + 
  facet_wrap(~chrom,nrow=1,scales="free_x") +
  ylab("SNVs / kb") +
  xlab("Physical position (Mb)") +
  #ggtitle("KD vs. KD")+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,100),expand = c(0.01,0))


divp <- cowplot::plot_grid(p1,p2,p3,nrow=3,align = "v",axis = "lr",labels = c("a","b","c"))

ggsave(divp,filename = "../../figures/SF13_rg_divergence.png",width = 7,height = 5.5,device = "png",units = "in",bg ="white",dpi = 900)
