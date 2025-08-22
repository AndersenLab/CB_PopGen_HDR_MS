rm(list=ls())

library(readxl)
library(geosphere)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)
library(reshape2)





# Load input data. File was downloaded on April 2nd 2024.
Cb_raw_csv<-read.csv("../../data/20250626_c_briggsae_strain_data.csv") 
# # Or load the data from google drive 
# raw_xlsx <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/edit")



plot_isotype_dist <- function(meta_sheet, condordance_groups, title_text,
                              bins = 20) {
  
  df <- meta_sheet %>%
    filter(strain %in% condordance_groups$strain) %>%
    select(-isotype) %>%
    left_join(condordance_groups, by = "strain")
  
  coords <- df %>%
    select(strain, longitude, latitude) %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- geosphere::distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- reshape2::melt(dist_mat, varnames = c("strain1", "strain2"), value.name = "geo_m") %>%
    filter(!is.na(geo_m)) %>%
    mutate(geo_distance = geo_m / 1000) %>%
    select(strain1, strain2, geo_distance)
  
  isotype_pairs <- df %>%
    select(strain, isotype) %>%
    filter(!is.na(isotype)) %>%
    group_by(isotype) %>%
    filter(n() > 1) %>%
    summarise(pairs = list(apply(combn(strain, 2), 2, paste, collapse = "_"))) %>%
    unnest(pairs) %>%
    pull(pairs)
  
  all_dist <- all_dist %>%
    mutate(pair1 = paste(strain1, strain2, sep = "_"),
           pair2 = paste(strain2, strain1, sep = "_"))
  dist_iso <- all_dist %>%
    filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>%
    pull(geo_distance)
  
  p <- ggplot(data.frame(dist = dist_iso), aes(x = dist)) +
    geom_histogram(bins = bins, fill = "#53886C", color = "black") +
    stat_bin(bins = bins, geom = "text",
             aes(label = after_stat(count)), vjust = -0.5, size = 3) +
    labs(title = title_text,
         x = "Distance (km)",
         y = "Frequency") +
    theme_classic()
  
  
  dist_out <- all_dist %>%
    filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>% 
    select(strain1, strain2, geo_distance) %>% 
    arrange(desc(geo_distance))
  
  
  return(list(plot=p,
              strain_list=dist_out))
  
}





# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE)

Cb_list <- plot_isotype_dist(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = expression(paste(italic("C. briggsae"), " - Haversine distance between strains within an isotype"))
)

Cb_list$plot
# ggsave("Cb_Haversine_distance_between_strains_main.pdf", plot = Cb_list$plot, width = 7.5, height = 4, units = "in")
ggsave("Cb_Haversine_distance_between_strains_main.pdf", plot = Cb_list$plot, width = 7, height = 4, units = "in")



### Add isotype info
Cb_condordance_groups_tmp<-Cb_condordance_groups_raw %>% 
  select(strain,isotype) %>% 
  rename(strain1=strain)

Cb_geo_distance_within_same_isotype<- Cb_list$strain_list %>% 
  left_join(Cb_condordance_groups_tmp, by =c("strain1") )
  
write.csv(Cb_geo_distance_within_same_isotype,
          "../../processed_data/Geo_info/Cb_geo_distance_within_same_isotype.csv",
          quote = FALSE,
          row.names = FALSE) 




######### How many isotypes have strains greater than 100 km from each other? 
######### Defined as global isotypes ########
Cb_geo_distance_within_same_isotype<-read.csv("../../processed_data/Geo_info/Cb_geo_distance_within_same_isotype.csv")

Cb_df_geo_over100km<-Cb_geo_distance_within_same_isotype %>% 
  filter(geo_distance >100)

Cb_100km_isotype<-unique(Cb_df_geo_over100km$isotype)
n_Cb_100km_isotype<-length(Cb_100km_isotype)
n_Cb_100km_isotype







######## add geo of isotypes > 100 km
indep_info_geo_for_each_strain<-read.csv("../../processed_data/Geo_info/indep_info_geo_for_each_strain.csv")

Cb_over_100km_df<-indep_info_geo_for_each_strain %>% 
  filter(isotype %in% Cb_100km_isotype) %>% 
  na.omit()
length(Cb_over_100km_df$isotype %>% unique())
## 38





#### Setting 100 km as the threshold for defining "cosmopolitan" samples


Cb_df_geo_over100km<-Cb_geo_distance_within_same_isotype %>% 
  filter(geo_distance >100)

length(unique(Cb_df_geo_over100km$isotype))
#38

Cb_cosmopolitan_isotype<-unique(Cb_df_geo_over100km$isotype)
Cb_cosmopolitan_isotype

write.table(Cb_cosmopolitan_isotype,
            "../../processed_data/Geo_info/Cb_cosmopolitan_isotype.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE) 










###################################
####### plot_figure_S1-c-d #########
##########################################

library(readxl)
library(geosphere)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)
library(reshape2)



plot_fs1 <- function(meta_sheet, condordance_groups, title_text,
                     bins = 22) {
  library(dplyr)
  library(ggplot2)
  
  df <- meta_sheet %>%
    filter(strain %in% condordance_groups$strain) %>%
    select(-isotype) %>%
    left_join(condordance_groups, by = "strain")
  
  coords <- df %>%
    select(strain, longitude, latitude) %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- geosphere::distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- reshape2::melt(dist_mat,
                             varnames = c("strain1", "strain2"),
                             value.name = "geo_m") %>%
    filter(!is.na(geo_m)) %>%
    mutate(geo_distance = geo_m / 1000) %>%
    select(strain1, strain2, geo_distance)
  
  isotype_pairs <- df %>%
    select(strain, isotype) %>%
    filter(!is.na(isotype)) %>%
    group_by(isotype) %>%
    filter(n() > 1) %>%
    summarise(pairs = list(apply(combn(strain, 2), 2, paste, collapse = "_"))) %>%
    unnest(pairs) %>%
    pull(pairs)
  
  all_dist <- all_dist %>%
    mutate(pair1 = paste(strain1, strain2, sep = "_"),
           pair2 = paste(strain2, strain1, sep = "_"))
  
  dist_iso <- all_dist %>%
    filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>%
    pull(geo_distance)
  
  # —— 新：手动分箱 —— 
  # 1. 中间 [100,5000] 区间等距分成 bins-1 段
  mid_breaks <- seq(100, 5000, length.out = bins - 1)
  # 2. 完整的 break 向量，加上 -Inf 和 Inf
  breaks <- c(-Inf, mid_breaks, Inf)
  # 3. 生成对应的标签
  labels <- c(
    "<100",
    paste0(round(mid_breaks[-length(mid_breaks)]), "–", round(mid_breaks[-1])),
    ">5000"
  )
  # 4. 用 cut() 生成一个有序因子
  binned <- cut(dist_iso,
                breaks = breaks,
                labels = labels,
                include.lowest = TRUE,
                right = FALSE)
  
  df_bin <- data.frame(bin = binned) %>%
    count(bin) %>%
    # 保证所有 bin 都在表里，即使 count 为 0
    complete(bin = factor(labels, levels = labels), fill = list(n = 0))
  
  # —— 用 geom_col() 绘制条形图 —— 
  p <- ggplot(df_bin, aes(x = bin, y = n)) +
    geom_col(fill = "grey30", color = NA) +
    geom_text(aes(label = n), vjust = -0.5, size = 1.8) +
    labs(title = title_text,
         x = "Distance (km)",
         y = "Frequency") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
          axis.text.y = element_text(size = 6))+
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) 
  
  
  dist_out <- all_dist %>%
    filter(pair1 %in% isotype_pairs | pair2 %in% isotype_pairs) %>%
    select(strain1, strain2, geo_distance) %>%
    arrange(desc(geo_distance))
  

  
  # 在脚本开头安装并加载 ggbreak
  # install.packages("ggbreak")
  library(ggbreak)
  
  # ...（前面的代码不变，生成 p）...
  
  # 对 p 进行断轴处理：断开 7000 到 14000 的区间
  p2 <- p +
    scale_y_break(c(6200, 14500),                 # 设置断轴区间
                  expand = c(0, 100))+
    scale_y_continuous(limits = c(0, 16700),
                       expand = c(0, 0))+
    theme(
      # 把右侧的刻度、文字和标题全都去掉
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.line.y.right  = element_blank()
    )+
    theme(
      plot.title      = element_text(size = 8, hjust = 0.5),  # 主标题
      axis.title.x    = element_text(size = 6,
                                     margin = margin(t = -6)),   # x轴标题
      axis.title.y    = element_text(size = 6,
                                     margin = margin(r = -6)),    # y轴标题
    )

  
  return(list(plot = p,
              plot_break = p2,
              strain_list = dist_out))
  
  
}






# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE)

Cb_fs1 <- plot_fs1(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = expression(paste(
    # italic("C. briggsae"), 
    "Pairwise Haversine distance between strains within an isotype"))
)

Cb_fs1$plot
Cb_fs1$plot_break

View(Cb_fs1$strain_list)

ggsave("FS1b_plot.png", 
       plot = Cb_fs1$plot, 
       width = 7,
       height = 3, 
       units = "in")

ggsave("FS1b_plot_break.png", 
       plot = Cb_fs1$plot_break, 
       width = 7,
       height = 3, 
       units = "in")





















#### Max Haversine distance between strains within the same isotype
library(dplyr)
library(ggplot2)
library(geosphere)
library(reshape2)
library(ggbreak)

plot_max_fs1 <- function(meta_sheet, condordance_groups, title_text, bins = 22) {
  # —— 原始数据预处理，与 plot_fs1 一样 ——  
  df <- meta_sheet %>%
    filter(strain %in% condordance_groups$strain) %>%
    select(-isotype) %>%
    left_join(condordance_groups, by = "strain") %>%
    filter(!is.na(isotype))
  
  coords <- df %>%
    select(strain, longitude, latitude) %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  # 计算全体 pairwise 距离矩阵
  dist_mat <- distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- melt(dist_mat,
                   varnames = c("strain1", "strain2"),
                   value.name = "geo_m") %>%
    filter(!is.na(geo_m)) %>%
    mutate(geo_distance = geo_m / 1000) %>%
    select(strain1, strain2, geo_distance)
  
  # —— 关键：只保留同一 isotype 的 pair，然后按 isotype 取最大距离 ——  
  # 1) 把 isotype 信息 join 到所有 pair 上
  iso_info <- df %>% select(strain, isotype)
  dist_iso <- all_dist %>%
    left_join(iso_info, by = c("strain1" = "strain")) %>%
    rename(isotype1 = isotype) %>%
    left_join(iso_info, by = c("strain2" = "strain")) %>%
    rename(isotype2 = isotype) %>%
    filter(!is.na(isotype1), !is.na(isotype2), isotype1 == isotype2) %>%
    # 2) 按 isotype 分组并计算最大距离
    group_by(isotype = isotype1) %>%
    summarise(max_distance = max(geo_distance)) %>%
    ungroup()
  
  # 提取所有 isotype 的最大距离向量
  dist_max <- dist_iso$max_distance
  
  # —— 分箱逻辑同 plot_fs1 ——  
  mid_breaks <- seq(100, 5000, length.out = bins - 1)
  breaks <- c(-Inf, mid_breaks, Inf)
  labels <- c(
    "<100",
    paste0(round(mid_breaks[-length(mid_breaks)]), "–", round(mid_breaks[-1])),
    ">5000"
  )
  binned <- cut(dist_max,
                breaks = breaks,
                labels = labels,
                include.lowest = TRUE,
                right = FALSE)
  
  df_bin <- tibble(bin = binned) %>%
    count(bin) %>%
    complete(bin = factor(labels, levels = labels), fill = list(n = 0))
  
  # —— 绘图与断轴——
  p <- ggplot(df_bin, aes(x = bin, y = n)) +
    geom_col(fill = "#53886C", color = NA) +
    geom_text(aes(label = n), vjust = -0.5, size = 1.8) +
    labs(title = title_text,
         x = "Max distance (km)",
         y = "Number of isotypes") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))+
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) 
  
  
  p2 <- p +
    scale_y_break(c(18,215),
                  expand = c(0, 1)) +
    scale_y_continuous(limits = c(0, 220), expand = c(0, 0)) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.line.y.right  = element_blank()
    )+
    theme(
      plot.title      = element_text(size = 8, hjust = 0.5),  # 主标题
      axis.title.x    = element_text(size = 6,
                                     margin = margin(t = -6)),   # x轴标题
      axis.title.y    = element_text(size = 6,
                                     margin = margin(r = -6)),    # y轴标题
    )
  
  return(list(
    plot     = p,
    plot_break = p2,
    max_distances = dist_iso
  ))
}





# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE)

Cb_fs1_max <- plot_max_fs1(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = expression(paste(
    # italic("C. briggsae"), 
    "Max Haversine distance between strains within an isotype"))
)

Cb_fs1_max$plot
Cb_fs1_max$plot_break


ggsave("FS1b_max_plot.png", 
       plot = Cb_fs1_max$plot, 
       width = 7,
       height = 3, 
       units = "in")

ggsave("FS1b_max_plot_break.png", 
       plot = Cb_fs1_max$plot_break, 
       width = 7,
       height = 3, 
       units = "in")

















#### With annotation Max Haversine distance between strains within the same isotype
plot_max_fs1_annot <- function(meta_sheet, condordance_groups, title_text) {
  # —— 数据预处理 ——  
  df <- meta_sheet %>%
    filter(strain %in% condordance_groups$strain) %>%
    select(-isotype) %>%
    left_join(condordance_groups, by = "strain") %>%
    filter(!is.na(isotype))
  
  coords <- df %>%
    select(strain, longitude, latitude) %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    mutate(longitude = as.numeric(longitude),
           latitude  = as.numeric(latitude)) %>%
    tibble::column_to_rownames("strain")
  
  # pairwise 距离矩阵
  dist_mat <- distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- melt(dist_mat,
                   varnames = c("strain1", "strain2"),
                   value.name = "geo_m") %>%
    filter(!is.na(geo_m)) %>%
    mutate(geo_distance = geo_m / 1000) %>%
    select(strain1, strain2, geo_distance)
  
  # 只保留同 isotype，计算每个 isotype 最大距离
  iso_info <- df %>% select(strain, isotype)
  dist_iso <- all_dist %>%
    left_join(iso_info, by = c("strain1" = "strain")) %>%
    rename(isotype1 = isotype) %>%
    left_join(iso_info, by = c("strain2" = "strain")) %>%
    rename(isotype2 = isotype) %>%
    filter(isotype1 == isotype2) %>%
    group_by(isotype = isotype1) %>%
    summarise(max_distance = max(geo_distance)) %>%
    ungroup()
  
  # —— 分箱并标记 cosmopolitan/Non-cosmopolitan ——  
  breaks <- c(-Inf, 1, 10, 50, 100, 500, 2000, Inf)
  labels <- c("< 1", "1-10", "10-50", "50-100", "100-500", "500-2000", "> 2000")
  
  df_bin <- dist_iso %>%
    mutate(
      bin = cut(max_distance,
                breaks = breaks,
                labels = labels,
                include.lowest = TRUE,
                right = FALSE),
      category = ifelse(max_distance > 100, "cosmopolitan", "Non-cosmopolitan")
    ) %>%
    count(bin, category) %>%
    complete(
      bin = factor(labels, levels = labels),
      category = c("cosmopolitan", "Non-cosmopolitan"),
      fill = list(n = 0)
    )
  
  # —— 计算每个 bin 的总数并绘图 ——  
  total_bin <- df_bin %>%
    group_by(bin) %>%
    summarise(total = sum(n))
  
  p <- ggplot(df_bin, aes(x = bin, y = n, fill = category)) +
    geom_col() +
    geom_text(
      data = total_bin,
      aes(x = bin, y = total, label = total),
      vjust = -0.3, size = 2.5, inherit.aes = FALSE
    ) +
    scale_fill_manual(
      breaks = c("Non-cosmopolitan", "cosmopolitan"),
      values = c("cosmopolitan" = "#53886C", "Non-cosmopolitan" = "grey"),
      guide = guide_legend(title = NULL)
    ) +
    labs(
      title = title_text,
      x = "Max distance (km)",
      y = "Number of isotypes"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      # legend.position = "inside",
      legend.position = c(0.8, 0.8),
      
      legend.key.size = unit(0.8, "lines"),
      legend.text    = element_text(size = 5),
      
      plot.title      = element_text(size = 8, hjust = 0.5),
      axis.title.x    = element_text(size = 6),
      axis.title.y    = element_text(size = 6)
    ) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8)))+
    scale_y_continuous(limits = c(0, 190), expand = c(0, 0))
    
    
    
  
  p2 <- p +
    scale_y_break(c(35, 170), expand = c(0, 1)) +
    scale_y_continuous(limits = c(0, 190), expand = c(0, 0)) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      legend.position = "none",
      plot.title      = element_text(size = 8, hjust = 0.5),
      axis.title.x    = element_text(size = 6, margin = margin(t = -6)),
      axis.title.y    = element_text(size = 6, margin = margin(r = -6))
    )
  
  return(list(
    plot         = p,
    plot_break   = p2,
    max_distances = dist_iso
  ))
}



# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE)

Cb_fs1_max_annot <- plot_max_fs1_annot(
  meta_sheet = Cb_raw_csv,
  condordance_groups = Cb_condordance_groups_raw,
  title_text = expression(paste(
    # italic("C. briggsae"), 
    "Max Haversine distance between strains within an isotype"))
)

Cb_fs1_max_annot$plot
Cb_fs1_max_annot$plot_break


ggsave("FS1b_max_plot_annot.png", 
       plot = Cb_fs1_max_annot$plot, 
       width = 7,
       height = 3, 
       units = "in",
       dpi = 600)

ggsave("FS1b_max_plot_annot.pdf", 
       plot = Cb_fs1_max_annot$plot, 
       width = 7,
       height = 3, 
       units = "in")

ggsave("FS1b_max_plot_break_annot.png", 
       plot = Cb_fs1_max_annot$plot_break, 
       width = 7,
       height = 3, 
       units = "in",
       dpi = 600)















############## strains counts ##########

# 
# meta_sheet = Cb_raw_csv
# concordance_groups = Cb_condordance_groups_raw
# title_text = expression(paste(
#   # italic("C. briggsae"), 
#   "Max Haversine distance between strains within an isotype"))

plot_strain_counts_per_bin <- function(meta_sheet, concordance_groups, title_text) {
  # —— 数据预处理：直接 join 出 isotype ——  
  df <- meta_sheet %>%
    inner_join(
      concordance_groups %>% select(strain, isotype),
      by = c("strain","isotype")
    )
  
  # —— 计算每个 isotype 的最大地理距离 ——  
  coords <- df %>%
    select(strain, longitude, latitude) %>%
    filter(!is.na(longitude), !is.na(latitude)) %>%
    mutate(
      longitude = as.numeric(longitude),
      latitude  = as.numeric(latitude)
    ) %>%
    tibble::column_to_rownames("strain")
  
  dist_mat <- distm(coords, fun = distHaversine)
  dist_mat[upper.tri(dist_mat, diag = TRUE)] <- NA
  rownames(dist_mat) <- colnames(dist_mat) <- rownames(coords)
  
  all_dist <- melt(dist_mat,
                   varnames = c("strain1", "strain2"),
                   value.name = "geo_m") %>%
    filter(!is.na(geo_m)) %>%
    mutate(geo_distance = geo_m / 1000) %>%
    select(strain1, strain2, geo_distance)
  
  iso_info <- df %>% select(strain, isotype)
  dist_iso <- all_dist %>%
    left_join(iso_info, by = c("strain1" = "strain")) %>%
    rename(isotype1 = isotype) %>%
    left_join(iso_info, by = c("strain2" = "strain")) %>%
    rename(isotype2 = isotype) %>%
    filter(isotype1 == isotype2) %>%
    group_by(isotype = isotype1) %>%
    summarise(max_distance = max(geo_distance)) %>%
    ungroup()
  
  # —— 分箱规则 ——  
  breaks <- c(-Inf, 1, 10, 50, 100, 500, 2000, Inf)
  labels <- c("< 1", "1-10", "10-50", "50-100", "100-500", "500-2000", "> 2000")
  
  # —— 给每个 strain 打 bin ——  
  strain_bins <- df %>%
    select(strain, isotype) %>%
    left_join(dist_iso, by = "isotype") %>%
    mutate(bin = cut(max_distance,
                     breaks = breaks,
                     labels = labels,
                     include.lowest = TRUE,
                     right = FALSE)) %>% 
    na.omit() %>% 
    mutate(category = ifelse(max_distance > 100, "cosmopolitan", "Non-cosmopolitan")
           )
  
  # —— 统计每个 bin 的 strain 数量 ——  
  strain_count <- strain_bins %>%
    count(bin, name = "n_strains") %>%
    complete(bin = factor(labels, levels = labels),
             fill = list(n_strains = 0)) %>% 
    mutate(category = ifelse(bin %in% c("100-500", "500-2000", "> 2000"),
                             "cosmopolitan", "Non-cosmopolitan")) %>% 
    mutate(category = factor(category,
                             levels = c("Non-cosmopolitan","cosmopolitan")))
  

  
    p_strains <- ggplot(strain_count, aes(x = bin, 
                                          y = n_strains,
                                          fill = category)) +
    geom_col() +
    geom_text(aes(label = n_strains),
              vjust = -0.3, size = 2.5) +
      scale_fill_manual(
        values = c(
          "Non-cosmopolitan" = "grey80",
          "cosmopolitan"     = "#53886C"
        ),guide = guide_legend(title = NULL)
      )+
    labs(
      title = title_text,
      x     = "Max distance (km)",
      y     = "Number of strains"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title  = element_text(size = 8, hjust = 0.5),
      axis.title  = element_text(size = 6),
      legend.position = c(0.8, 0.8),
      legend.key.size = unit(0.8, "lines"),
      legend.text    = element_text(size = 5)
    ) +
    scale_x_discrete(expand = expansion(add = c(0.8, 0.8))) +
    scale_y_continuous(limits = c(0, 700),
                       expand = c(0, 0))
  
  return(p_strains)
}






# add latest isotyoe info 
Cb_condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                                      header = TRUE)

Cb_fs1_max_annot_strains <- plot_strain_counts_per_bin(
  meta_sheet = Cb_raw_csv,
  concordance_groups = Cb_condordance_groups_raw,
  title_text = expression(paste(
    # italic("C. briggsae"), 
    "Max Haversine distance between strains within an isotype"))
)

Cb_fs1_max_annot_strains
















######### combine Cb_fs1_max_annot$plot_break and Cb_plot_UniStr$plot
library(patchwork)
p1 <- Cb_fs1_max_annot$plot
p2 <- Cb_fs1_max_annot_strains

combined <- ggpubr::ggarrange(ggpubr::ggarrange(p1,
                                                p2,
                                         ncol = 2, 
                                         labels = c("c","d"),
                                         widths = c(0.5, 0.5)))
combined
saveRDS(combined,"../../processed_data/assemble_figure_S1/geo_distance_plots.RData")


ggsave("raw_FS1b_combined.png",
       plot = combined,
       width = 7,
       height = 2.5,
       units = "in",
       dpi = 600)
ggsave("raw_FS1b_combined.pdf",
       plot = combined,
       width = 7,
       height = 2.5,
       units = "in")









