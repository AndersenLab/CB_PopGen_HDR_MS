library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(dendextend)
library(ComplexHeatmap)
library(ape) 
library(circlize)


hdrs <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_allStrain_5kbclust_20250930.tsv") 
hdrs_unt_rg <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_otherRG_UNT_5kbclust_20250930.tsv") 
bins <- readr::read_tsv("../../processed_data/HDRs/QX1410_genomic_windows.1kb.bed",col_names = c("CHROM","binStart","binEnd")) 
variants_QX <- readr::read_tsv("../../processed_data/HDRs/Tropical.variant_counts.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN")) %>% dplyr::mutate(source="QX1410") %>% dplyr::filter(!STRAIN=="QX1410")
variants_rest <- readr::read_tsv("../../processed_data/HDRs/Other_RG.variant_counts.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN"))#%>% dplyr::filter(STRAIN!="QX1410" & STRAIN!="JU2536")
variants_nrg <- readr::read_tsv("../../processed_data/HDRs/CB_all_strain_vc.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN"))%>% dplyr::filter(STRAIN!="QX1410" & STRAIN!="JU2536")
lineages <- readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRGeage_GeoLocAdmCol_20250909.tsv") %>%
  dplyr::mutate(sublineage_color=ifelse(Sublineage=="TC","#ff0000",sublineage_color)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="TC","TT",Sublineage))  %>% 
  dplyr::mutate(REF=ifelse(Lineage=="Tropical","QX1410",
                           ifelse(Lineage=="Temperate","JU2536",
                                  ifelse(Lineage=="TH","NIC1660",
                                         ifelse(Lineage=="AD","ECA2670",
                                                ifelse(Lineage=="KD","JU1348",
                                                       ifelse(Lineage=="TD1","BRC20530",
                                                              ifelse(Lineage=="TD2","BRC20492",NA))))))))
variants_rest_byGroup <- variants_rest %>% 
  dplyr::left_join(lineages %>% dplyr::select(isotype,REF),by=c("STRAIN"="isotype")) %>% 
  dplyr::rename(source=REF) 
all_variants <- rbind(variants_QX,variants_rest_byGroup)
hdrs_unt_all <- rbind(hdrs %>% dplyr::filter(source=="QX1410"),hdrs_unt_rg)
variants_nrg_byGroup <- variants_nrg %>% dplyr::left_join(hdrs_unt_all %>% dplyr::select(STRAIN,source) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>% dplyr::filter(!is.na(source))


dt_variants <- as.data.table(all_variants)
dt_hdrs <- as.data.table(hdrs_unt_all)

setnames(dt_variants, old = c("START_BIN", "END_BIN"), new = c("start", "end"))
setnames(dt_hdrs, old = c("minStart", "maxEnd"), new = c("start", "end"))

setkey(dt_variants, source, STRAIN, CHROM, start, end)
setkey(dt_hdrs, source, STRAIN, CHROM,  start, end)

overlap_result <- foverlaps(dt_variants, dt_hdrs, type = "any", nomatch = 0L)

dt_variants[, in_hd := FALSE]
overlap_bins <- overlap_result[, .(source, STRAIN, CHROM, start = i.start, end = i.end)]
dt_variants[overlap_bins, in_hd := TRUE, on = .(source, STRAIN, CHROM, start, end)]

countSummary_nr <- as.data.frame(dt_variants) %>%
  dplyr::group_by(STRAIN,in_hd) %>%
  dplyr::summarise(tot_var=sum(COUNT),source=dplyr::first(source))

variantSummary_nr <- countSummary_nr %>%
  tidyr::pivot_wider(names_from = in_hd, values_from = tot_var, names_prefix = "tot_var_") %>%
  dplyr::rename(non_divergent_variants=tot_var_FALSE,divergent_variants=tot_var_TRUE) %>%
  dplyr::mutate(genome_wide_variants=non_divergent_variants+divergent_variants,divergent_variants_prop=divergent_variants/genome_wide_variants) 

covSpans_nr <- as.data.frame(dt_variants) %>%
  dplyr::group_by(STRAIN,in_hd) %>%
  dplyr::summarise(span=n()*1e3)

spanSummary_nr <- covSpans_nr %>%
  tidyr::pivot_wider(names_from = in_hd, values_from = span, names_prefix = "tot_bases_") %>%
  dplyr::rename(non_divergent_span=tot_bases_FALSE,divergent_span=tot_bases_TRUE) %>%
  dplyr::mutate(genome_span=non_divergent_span+divergent_span,divergent_span_prop=divergent_span/genome_span) 

allSummary_nr <- variantSummary_nr %>%
  dplyr::left_join(spanSummary_nr,by="STRAIN") %>%
  dplyr::filter(source!=STRAIN)

meanRGSummary_nr <- allSummary_nr %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(mean_var_prop = mean(divergent_variants_prop),
                   mean_var_prop_sd=sd(divergent_variants_prop),
                   min_var_prop=min(divergent_variants_prop),
                   max_var_prop=max(divergent_variants_prop),
                   mean_span_prop=mean(divergent_span_prop),
                   mean_span_prop_sd=sd(divergent_span_prop),
                   min_span_prop=min(divergent_span_prop),
                   max_span_prop=max(divergent_span_prop)) %>%
  dplyr::ungroup()

custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E6D700",
  "#A65628", "#F781BF", "#444444", "#90EE90", "#00CED1", "#8DA0CB")
lineage_cols <-lineages %>% dplyr::select(Lineage, lineage_color) %>%
  dplyr::distinct(Lineage, lineage_color) %>%
  tibble::deframe()

leg_source <- ggplot() +
  geom_point(data=allSummary_nr  %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
             aes(x = divergent_span_prop * 100, y = divergent_variants_prop * 100, color = Lineage,shape="RG")) +
  geom_point(data=allSummary_nr  %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
             aes(x = divergent_span_prop * 100, y = divergent_variants_prop * 100, color = Lineage,shape="QX1410")) +
  scale_color_manual(values=lineage_cols) + theme_bw() +
  labs(color = "Relatedness\ngroup",shape = "Reference\ngenome")+
  scale_shape_manual(values = c(0,1)) 



p1_nr <- ggplot() +
  geom_point(data=allSummary_nr %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
             aes(x = divergent_span_prop * 100, y = divergent_variants_prop * 100, color = Lineage),shape=1,stroke=0.8) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none',
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    plot.margin = margin(5.5, 1, 1, 1)
  ) +
  scale_color_manual(values=lineage_cols) +
  scale_y_continuous(breaks = seq(0, 80, 10),limits = c(0,85)) +
  scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,24)) +
  labs(color = "Relatedness group")


p2_nr <- ggplot() +
  geom_point(data=meanRGSummary_nr %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
             aes(x = mean_span_prop * 100, y = mean_var_prop * 100, color = Lineage)) +
  geom_errorbar(data=meanRGSummary_nr %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
                aes(x = mean_span_prop * 100,
                   ymin = (mean_var_prop - mean_var_prop_sd) * 100,
                   ymax = (mean_var_prop + mean_var_prop_sd) * 100,
                   color = Lineage),
               width = 0) + 
  geom_errorbarh(data=meanRGSummary_nr %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
                 aes(y = mean_var_prop * 100,
                     xmin = (mean_span_prop - mean_span_prop_sd) * 100,
                     xmax = (mean_span_prop + mean_span_prop_sd) * 100,
                     color = Lineage),
                 height = 0) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    plot.margin = margin(5.5, 1, 1, 1),
    legend.position = 'none') +
  scale_color_manual(values=lineage_cols) +
  scale_y_continuous(breaks = seq(0, 80, 10),limits = c(0,85)) +
  scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,24)) +
  labs(color = "Relatedness group")

legend <- get_legend(leg_source)

dt_variants <- as.data.table(variants_nrg_byGroup)
dt_hdrs <- as.data.table(hdrs)

setnames(dt_variants, old = c("START_BIN", "END_BIN"), new = c("start", "end"))
setnames(dt_hdrs, old = c("minStart", "maxEnd"), new = c("start", "end"))

setkey(dt_variants, source, STRAIN, CHROM, start, end)
setkey(dt_hdrs, source, STRAIN, CHROM,  start, end)

overlap_result <- foverlaps(dt_variants, dt_hdrs, type = "any", nomatch = 0L)

dt_variants[, in_hd := FALSE]
overlap_bins <- overlap_result[, .(source, STRAIN, CHROM, start = i.start, end = i.end)]
dt_variants[overlap_bins, in_hd := TRUE, on = .(source, STRAIN, CHROM, start, end)]

countSummary <- as.data.frame(dt_variants) %>%
  dplyr::group_by(STRAIN,in_hd) %>%
  dplyr::summarise(tot_var=sum(COUNT),source=dplyr::first(source))

variantSummary <- countSummary %>%
  tidyr::pivot_wider(names_from = in_hd, values_from = tot_var, names_prefix = "tot_var_") %>%
  dplyr::rename(non_divergent_variants=tot_var_FALSE,divergent_variants=tot_var_TRUE) %>%
  dplyr::mutate(genome_wide_variants=non_divergent_variants+divergent_variants,divergent_variants_prop=divergent_variants/genome_wide_variants) 

covSpans <- as.data.frame(dt_variants) %>%
  dplyr::group_by(STRAIN,in_hd) %>%
  dplyr::summarise(span=n()*1e3)

spanSummary <- covSpans %>%
  tidyr::pivot_wider(names_from = in_hd, values_from = span, names_prefix = "tot_bases_") %>%
  dplyr::rename(non_divergent_span=tot_bases_FALSE,divergent_span=tot_bases_TRUE) %>%
  dplyr::mutate(genome_span=non_divergent_span+divergent_span,divergent_span_prop=divergent_span/genome_span) 

allSummary <- variantSummary %>%
  dplyr::left_join(spanSummary,by="STRAIN") %>%
  dplyr::filter(!is.na(divergent_variants)) %>%
  dplyr::filter(source!=STRAIN)

meanRGSummary <- allSummary %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(mean_var_prop = mean(divergent_variants_prop),
                   mean_var_prop_sd=sd(divergent_variants_prop),
                   min_var_prop=min(divergent_variants_prop),
                   max_var_prop=max(divergent_variants_prop),
                   mean_span_prop=mean(divergent_span_prop),
                   mean_span_prop_sd=sd(divergent_span_prop),
                   min_span_prop=min(divergent_span_prop),
                   max_span_prop=max(divergent_span_prop)) %>%
  dplyr::ungroup()

p1 <- ggplot(allSummary %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype"))) +
  geom_point(aes(x = divergent_span_prop * 100, y = divergent_variants_prop * 100, color = Lineage),shape=0,stroke=0.8) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin = margin(5.5, 1, 1, 1),
    legend.position = "none"
  ) +
  scale_color_manual(values=lineage_cols) +
  scale_y_continuous(breaks = seq(0, 80, 10),limits = c(0,85)) +
  scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,24)) 


p2 <- ggplot() +
  geom_point(data=meanRGSummary %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
             aes(x = mean_span_prop * 100, y = mean_var_prop * 100, color = Lineage),shape=15) +
  geom_errorbar(data=meanRGSummary %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
                aes(x = mean_span_prop * 100,
                    ymin = (mean_var_prop - mean_var_prop_sd) * 100,
                    ymax = (mean_var_prop + mean_var_prop_sd) * 100,
                    color = Lineage),
                width = 0) + 
  geom_errorbarh(data=meanRGSummary %>% dplyr::left_join(lineages %>% dplyr::select(isotype,Lineage,lineage_color),by=c("source"="isotype")),
                 aes(y = mean_var_prop * 100,
                     xmin = (mean_span_prop - mean_span_prop_sd) * 100,
                     xmax = (mean_span_prop + mean_span_prop_sd) * 100,
                     color = Lineage),
                 height = 0) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(5.5, 1, 1, 1),
    legend.position = "none") +
  scale_color_manual(values=lineage_cols) +
  scale_y_continuous(breaks = seq(0, 80, 10),limits = c(0,85)) +
  scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,24)) +
  labs(color = "Relatedness group")

main_plot <- plot_grid(p1, p1_nr,p2, p2_nr, nrow = 2, align = "v", axis = "lr", rel_widths = c(1, 1), labels = c("a","b","c","d"),label_y = 1.01) 
#main_plot_with_legend <- plot_grid(main_plot, legend, nrow = 1, rel_widths = c(1, 0.2))

padded_plot <- plot_grid(NULL, main_plot, NULL, ncol = 3, rel_widths = c(0.04, 1, 0.01))
padded_plot2 <- plot_grid(NULL, padded_plot, NULL, nrow = 3, rel_heights = c(0.03, 1, 0.03))


full_plot <- ggdraw(padded_plot2) +
  draw_label("Percent genome span of hyper-divergent regions", x = 0.5, y = 0.01, vjust = 0, size = 12) +
  draw_label("Percent of variants in hyper-divergent regions", x = 0.01, y = 0.5, angle = 90, vjust = 1, size = 12)

full_plot_wleg <- plot_grid(full_plot, legend, nrow = 1, rel_widths = c(1, 0.18))

ggsave(plot = full_plot_wleg, filename = "../../figures/EDF5_propVC_byIsotype_20250930.png",width = 7,height = 6,dpi = 600,device = 'png',bg = "white")

summary_stats_perGroup <- rbind(meanRGSummary %>% dplyr::mutate(relative_to="QX1410"),meanRGSummary_nr %>% dplyr::mutate(relative_to="Relatedness Group"))

write.table(summary_stats_perGroup,file = "../../supplementary_data/SD7_summaryStats_20250930.tsv",sep = "\t",quote = F,row.names = F)

admix_color <- data.frame(
  letter = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
             "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U","admixed"),
  color = c("#4B0401", "#DA030F", "#FF6C93", "#D45602", "#563901",
            "#FFE5CA", "#FFB914", "#FFDA90", "#77C002", "#05E51B",
            "#02491E", "#05DEA2", "#82FFFA", "#01A4B1", "#5CB7FF",
            "#000F2D", "#0340B9", "#9C87FF", "#CBB1FF", "#F479FF", "#5F0158","gray80"),
  stringsAsFactors = FALSE
)

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
                "Cosmopolitan" = "gray30",
                "unknown" = 'grey')

df_colors <- data.frame(unname(geo_colors),names(geo_colors)) %>% dplyr::rename(color=`unname.geo_colors.`, geo=`names.geo_colors.`)

admix <- readr::read_tsv(file="../../processed_data/genetic_similarity_and_admixutre/non_admixed_isotypes.txt") %>%
  dplyr::select(samples,cluster)

geo <- readr::read_csv(file="../../processed_data/genetic_similarity_and_admixutre/Cb_indep_isotype_info_geo.csv") %>%
  dplyr::left_join(df_colors,by=c("geo")) %>%
  dplyr::mutate(abslat=abs(lat)) %>%
  dplyr::left_join(admix,by=c("isotype"="samples")) %>%
  dplyr::mutate(cluster=ifelse(is.na(cluster),"admixed",cluster)) %>%
  dplyr::left_join(admix_color,by=c("cluster"="letter")) %>%
  dplyr::rename(subpop=cluster)


conc <- readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/gtcheck.tsv")

isos <- readr::read_tsv(file="../../processed_data/genetic_similarity_and_admixutre/isotype_groups.tsv") %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(count=n())

lineages <- readr::read_tsv("../../processed_data/genetic_similarity_and_admixutre/isotype_byRG_GeoLocAdmCol_20250909.tsv") %>%
  dplyr::mutate(sublineage_color=ifelse(Sublineage=="TC","#ff0000",sublineage_color)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="TC","TT",Sublineage)) 

annotation_df <-as.data.frame(geo %>% dplyr::select(geo,abslat,subpop))
annotation_df$abslat <- as.numeric(annotation_df$abslat)
rownames(annotation_df) <- geo$isotype

QX1410_group <- c(unique((hdrs %>% dplyr::filter(source=="QX1410"))$STRAIN),"QX1410")

concordance_matrix <- conc %>%
  dplyr::filter(i %in% isos$isotype & j %in% isos$isotype) %>%
  dplyr::mutate(concordance = (sites - discordance) / sites) %>%
  dplyr::select(i, j, concordance) %>%
  dplyr::bind_rows(
    conc %>%
      dplyr::filter(i %in% isos$isotype & j %in% isos$isotype) %>%
      dplyr::mutate(concordance = (sites - discordance) / sites) %>%
      dplyr::select(i = j, j = i, concordance)
  ) %>%
  tidyr::pivot_wider(names_from = i, values_from = concordance) %>%
  tibble::column_to_rownames("j") %>%
  as.matrix()
all_iso <- sort(unique(isos$isotype))
concordance_matrix <- concordance_matrix[all_iso, all_iso]
diag(concordance_matrix) <- 1

filtered_mat <- concordance_matrix[QX1410_group, QX1410_group]

abslat_vec <- annotation_df$abslat
names(abslat_vec) <- rownames(annotation_df)
abslat_vec <- abslat_vec[rownames(filtered_mat)]  # align order

geo_vec <- annotation_df$geo
names(geo_vec) <- rownames(annotation_df)
geo_vec <- geo_vec[rownames(filtered_mat)] 

adm_vec <- annotation_df$subpop
names(adm_vec) <- rownames(annotation_df)
adm_vec <- adm_vec[rownames(filtered_mat)] 

subpop_levels <- setdiff(unique(adm_vec), "admixed")  # all clusters except admixed
adm_vec <- factor(adm_vec, levels = c(sort(subpop_levels), "admixed"))

admix_col_vec <- admix_color$color
names(admix_col_vec) <- admix_color$letter

# Color function for abslat (continuous gradient)
lat_col_fun <- circlize::colorRamp2(
  seq(0, 60, length.out = 5),
  c("#D73027", "#EDC948" , "#FFFFE0", "#4575B4", "#000080"),
)

abslat_annot <- ComplexHeatmap::HeatmapAnnotation(
  abslat = abslat_vec,
  col = list(abslat = lat_col_fun),
  annotation_name_side = "left",
  annotation_legend_param = list(
    abslat = list(title = "Absolute\nlatitude")
  ),
  which = "column",
  height = grid::unit(1.2, "cm")
)

geo_annot <- ComplexHeatmap::rowAnnotation(
  geo = geo_vec,
  col = list(geo = geo_colors),
  annotation_name_side = "top",
  annotation_legend_param = list(
    geo = list(title = "Geographic\nregion")
  ),
  width = grid::unit(1.2, "cm")
)
phylo_vals <- as.vector(concordance_matrix)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]  # remove NAs

# Set the breakpoints manually
breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)

# Define the color mapping with white at the median
phylo_col_fun <- circlize::colorRamp2(
  breaks,
  #seq(min(phylo_vals, na.rm = TRUE), max(phylo_vals, na.rm = TRUE), length.out = 5),
  c("#4575B4", "#74ADD1", "#FFFFE0","#F4D166","#D73027")
)

p_heatmap <- ComplexHeatmap::Heatmap(
  filtered_mat,
  name = "Genetic\nsimilarity",
  col = phylo_col_fun,
  cluster_rows = T,
  cluster_columns = T,
  clustering_distance_rows = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  show_row_names = F,
  show_column_names = F,
  right_annotation = geo_annot,
  bottom_annotation = abslat_annot,
)

QX1410_stats <- allSummary %>%
  dplyr::filter(source=="QX1410") %>%
  dplyr::left_join(geo,by=c("STRAIN"="isotype")) %>%
  dplyr::left_join(rbind(conc %>% dplyr::filter(i=="QX1410") %>% dplyr::select(discordance,sites,isotype=j),
                         conc %>% dplyr::filter(j=="QX1410") %>% dplyr::select(discordance,sites,isotype=i)) %>%
                     dplyr::distinct(isotype,.keep_all = T), by=c("STRAIN"="isotype")) %>%
  dplyr::mutate(concordance=(sites-discordance)/sites) %>%
  dplyr::left_join(lineages %>% dplyr::select(isotype,Sublineage,sublineage_color),by=c("STRAIN"="isotype"))
  

vals_rescaled <- scales::rescale(breaks)

pqx1 <- ggplot(QX1410_stats) +
  geom_point(aes(x = divergent_span_prop * 100, y = divergent_variants_prop * 100, fill = concordance),shape=22,color="grey20",stroke=0.3,size=2) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.81, 0.35),
    axis.title=element_text(size=10)
  ) +
  # scale_y_continuous(breaks = seq(0, 60, 5),limits = c(0,62)) +
  # scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,12)) +
  scale_fill_gradientn(
    colors = c("#4575B4", "#74ADD1", "#FFFFBF", "#F4D166", "#D73027"),
    values = vals_rescaled,
    limits = range(0.8250383,1),
    oob = scales::squish,
    name = "Genetic\nsimilarity\nwith QX1410\n"
  )+
  xlab("Percent genome span of hyper-divergent regions") +
  ylab("Percent of variants in hyper-divergent regions")

breaks_lat <- c(0,11.75, 23.5, 35.25,47, 60)

QX1410_stats$geo <- factor(QX1410_stats$geo, levels = names(geo_colors))

pqx2 <- ggplot(QX1410_stats) +
  geom_point(
    aes(x = divergent_span_prop * 100,
        y = divergent_variants_prop * 100,
        fill = geo),
    shape = 22, color = "grey20", stroke = 0.3, size = 2
  ) +
  scale_fill_manual(values = geo_colors) +  # <- map to your named vector
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.81, 0.25),
    axis.title = element_text(size = 9.5)
  ) +
  xlab("Percent genome span of hyper-divergent regions") +
  ylab("Percent of variants in hyper-divergent regions")

p_heatmap_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(p_heatmap,  merge_legend = TRUE))
pqx <- plot_grid(pqx1,pqx2,ncol=2,rel_widths = c(1,1), labels = c("b","c"), align = "h", axis = "tb")
qx_heat <- plot_grid(p_heatmap_grob,pqx,ncol=1, labels = c("a",NA), align = "v", axis = "lr",rel_heights = c(1.4,1))

#pqx1 + coord_cartesian(xlim=c(5,NA),ylim=c(0,6.5))
outliers <- QX1410_stats %>% dplyr::filter(divergent_variants_prop <= 0.065 & divergent_span_prop >=0.05)

fold_average <- allSummary %>%
  dplyr::filter(source=="QX1410") %>%
  dplyr::mutate(divergent_variant_density=divergent_variants/divergent_span,genome_wide_variant_density=genome_wide_variants/genome_span) %>%
  dplyr::mutate(fold_change=divergent_variant_density/genome_wide_variant_density) %>%
  dplyr::filter(STRAIN=="QX1410")

size_summary <- hdrs %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarize(total_divSize = sum(divSize, na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(total_divSize)) %>%
  dplyr::mutate(prop=total_divSize/106184000)

meanSize <- mean(hdrs$divSize) /1e3
minSize <- min(hdrs$divSize) /1e3
maxSize <- max(hdrs$divSize) /1e3

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

collapsed <- plyr::ldply(getRegFreq(hdrs %>% dplyr::group_split(CHROM)),data.frame) %>% dplyr::mutate(divSize=maxEnd-minStart)
cov_wg <- sum(collapsed$divSize) / 106184000
wg_tot <- nrow(collapsed)
wg_mean <- mean(collapsed$divSize) / 1e3
wg_min <- min(collapsed$divSize) / 1e3
wg_max <- max(collapsed$divSize) / 1e3

collapsed_bysource <- plyr::ldply(getRegFreq(hdrs %>%dplyr::filter() %>% dplyr::group_split(CHROM,source)),data.frame) %>% dplyr::mutate(divSize=maxEnd-minStart)

bySource_summary <- collapsed_bysource %>%
  dplyr::group_by(source) %>%
  dplyr::summarize(
  tot_reg = n(),
  mean_divSize = mean(divSize, na.rm = TRUE) /1e3,
  meadian_divSize = median(divSize,na.rm =TRUE)/1e3,
  min_divSize = min(divSize, na.rm = TRUE)/1e3,
  max_divSize = max(divSize, na.rm = TRUE)/1e3,
  percent_cov = sum(divSize, na.rm = TRUE) / 106184000
)

strain_counts <- hdrs %>%
  dplyr::group_by(STRAIN, source) %>%
  dplyr::summarize(n_regions = dplyr::n(), perc_cov=sum(divSize)/106184000) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(`Relatedness\ngroup`=ifelse(source=="QX1410","Tropical","Other"))


# Plot
p_nreg <- ggplot(strain_counts,  aes(x = forcats::fct_reorder(STRAIN, n_regions), y = n_regions, fill = `Relatedness\ngroup`)) +
  geom_col() +
  theme_bw() +
  labs(x = "695 isotype reference strains",y = "Number of HDRs",  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'inside',
        legend.position.inside = c(0.1,0.75),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color="grey",linetype = "dashed")) +
  scale_y_continuous(expand = c(0,0))

p_cov <- ggplot(strain_counts,  aes(x = forcats::fct_reorder(STRAIN, perc_cov*100), y = perc_cov*100, fill = `Relatedness\ngroup`)) +
  geom_col() +
  theme_bw() +
  labs(x = "695 isotype reference strains",y = "Percent genome covered",  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color="grey",linetype = "dashed"),
        legend.position = 'none') +
  scale_y_continuous(expand = c(0,0))

p_cov_nreg <- cowplot::plot_grid(p_nreg,p_cov,nrow=2,ncol=1, align = "v",labels = c("a","b"))

propQX_nreg <- nrow(strain_counts %>% dplyr::filter(n_regions>=300 & source=="QX1410")) / nrow(strain_counts %>% dplyr::filter(n_regions>=300))
propQX_pcov <- nrow(strain_counts %>% dplyr::filter(perc_cov>=0.05 & source=="QX1410")) / nrow(strain_counts %>% dplyr::filter(perc_cov>=0.05))

ggsave(plot = p_cov_nreg, filename = "../../figures/SF10_nreg_byIsotype_20250930.png",width = 7.5,height = 5.5,dpi = 600,device = 'png')

bins_dt <- as.data.table(bins)
setnames(bins_dt, c("binStart", "binEnd"), c("start", "end"))
bins_dt[, id := .I]  # optional: keep track of bins

hdrs_dt <- as.data.table(hdrs)
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
  dplyr::mutate(freq=n_strains/695)

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

x_breaks_by_chrom <- bins_wFreq %>%
  group_by(CHROM) %>%
  summarise(x_breaks = list(seq(
    floor(min((start) / 1e6)),
    ceiling(max((end) / 1e6)),
    by = 5
  )))

p1 <- ggplot(hdrs_ordered) + 
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
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 5),expand = c(0, 0))

p2 <- ggplot(bins_wFreq) +
  geom_point(aes(x=(start+500)/1e6,y=freq),size=0.2,stroke = 0) +
  facet_wrap(~CHROM,scales = 'free_x',nrow=1) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        strip.text = element_blank())  +
  ylab("Frequency")+
  xlab("Physical position (Mb)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x+(500/1e6))), by = 5),expand = c(0, 0))

p21 <- cowplot::plot_grid(p1,p2,nrow=2,ncol=1,rel_heights = c(1,0.3),align = "v", axis = "l", labels = c("a","b"))

