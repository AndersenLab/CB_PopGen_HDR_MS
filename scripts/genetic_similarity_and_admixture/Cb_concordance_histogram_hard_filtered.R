rm(list = ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)

source("../utilities.R")

gtcheck_raw <- read.delim("../../processed_data/genetic_similarity_and_admixutre/strain_gtcheck.txt") 

gtcheck<-gtcheck_raw %>% 
  dplyr::filter(!(i %in% c("MY681","ECA1146","JU356","ECA1503"))) %>% 
  dplyr::filter(!(j %in% c("MY681","ECA1146","JU356","ECA1503"))) %>% 
  dplyr::mutate(concordance = 1-(discordance/sites)) %>% 
  dplyr::select(i,j,concordance)

head(gtcheck)


th <- 0.9997

p_full <- ggplot(gtcheck, aes(x = concordance)) +
  geom_histogram(binwidth = 0.001, fill = "grey60", color = "grey30") +
  geom_vline(xintercept = th, color = "blue", size = 1) +
  labs(x = "Genetic similarity", y = "Number of Comparisons") +
  theme_bw(base_size = 7) +
  theme(
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold")
  )

xleft <- 0.999
xright <- 1.0
zoom_bin <- 0.00001

vals_zoom <- gtcheck$concordance
vals_zoom <- vals_zoom[!is.na(vals_zoom) & vals_zoom >= xleft & vals_zoom <= xright]

p_zoom <- ggplot(filter(gtcheck, !is.na(concordance) & concordance >= xleft), aes(x = concordance)) +
  geom_histogram(binwidth = zoom_bin, fill = "grey60", color = "grey30", alpha = 0.9, na.rm = TRUE) +
  annotate("rect", xmin = th, xmax = xright, ymin = -Inf, ymax = Inf, alpha = 0.4, fill = "dodgerblue", color = "blue") +
  coord_cartesian(xlim = c(xleft, xright), expand = FALSE) +
  labs(x = "Genetic similarity", y = "Number of Comparisons") +
  annotate("text", x = th - 0.00025, y = 4000,
           label = "the same isotype", color = "blue", size = 4, fontface = "bold") +
  theme_bw(base_size = 7)+
  theme(
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold")
  )

final <- plot_grid(p_full, p_zoom, ncol = 2, 
                   rel_widths = c(1, 1),
                   labels = c("a", "b"))

ggsave("../../figures/SF14_concordance_histogram.pdf",final, width = 7, height = 3, units = "in")

