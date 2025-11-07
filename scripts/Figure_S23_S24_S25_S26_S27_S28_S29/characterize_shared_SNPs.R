library(dplyr)
library(tidyr)
library(vcfR)
library(pheatmap)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(grid)
library(stringr)
library(ggplot2)
library(ggh4x)
library(data.table)

#read reference genome (1kb) bins
bins <- readr::read_tsv("../processed_data/QX1410_genomic_windows.1kb.bed", col_names = c("CHROM","start","end"))

#read reference genome chromosomal domain breakpoints
domains_raw <- readr::read_csv("../processed_data/chromosome_domain_Cbriggsae.csv") %>%
  dplyr::rename(CHROM=chrom,start=left,end=right) %>%
  dplyr::mutate(start=start*1e3,end=end*1e3)

#get reference coding sequences
cds <- readr::read_tsv("../processed_data/c_briggsae.QX1410_20250929.csq.gff", col_names = c("chrom","source","type","start","end","score","strand","phase","attributes")) %>%
  dplyr::filter(type=="CDS") %>%
  dplyr::select(chrom,start,end,strand,attributes) %>%
  tidyr::separate(attributes,into=c("ID","Parent","seqname","biotype","locus"),sep = ";") %>%
  dplyr::mutate(tran=gsub("Parent=transcript:","",Parent),seqname=gsub("sequence_name=","",seqname)) %>%
  dplyr::select(-ID,-biotype,-locus,-Parent) 

#get a list of Tropical isotype names
Tropical_samples <- readr::read_tsv("/../processed_data/isotype_byLineage_GeoLocAdmCol_20250909.tsv") %>% dplyr::filter(Lineage=="Tropical" & isotype!="QX1410") %>% dplyr::pull(isotype)

#read shared SNPs between Tropical and each relatedness group
KDshared <- readr::read_tsv("../processed_data/tropical_joins/Tropical_vs_KD.tsv") %>% dplyr::select(all_of(c("POS",Tropical_samples)))
ADshared <- readr::read_tsv("../processed_data/Tropical_vs_AD.tsv") %>% dplyr::select(all_of(c("POS",Tropical_samples)))
TD1shared <- readr::read_tsv("../processed_data/Tropical_vs_TD1.tsv") %>% dplyr::select(all_of(c("POS",Tropical_samples)))
TEMPshared <- readr::read_tsv("../processed_data/Tropical_vs_Temperate.tsv") %>% dplyr::select(all_of(c("POS",Tropical_samples)))
THshared <- readr::read_tsv("../processed_data/Tropical_vs_TH.tsv") %>% dplyr::select(all_of(c("POS",Tropical_samples)))

#read hyper-divergent region boundaries
hdrs <- readr::read_tsv("../processed_data/HDR_CB_allStrain_5kbclust_20250930.tsv") 

#read all tropical SNP sites, filter sites with >10% missing genotypes and <1% MAF
Trop_sites <- readr::read_tsv("../processed_data/Tropical_boolean_matrix_wHead.tsv") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::mutate(n_ones = rowSums(select(., -POS) == 1, na.rm = TRUE)) %>%
  dplyr::mutate(n_zeros = rowSums(select(., -POS) == 0, na.rm = TRUE)) %>%
  dplyr::filter(n_NA < 50) %>%
  dplyr::filter(n_zeros >= 5 & n_ones >= 5)

#read all Kerala Divergent SNP sites, filter sites with >10% missing genotypes
KD_sites <- readr::read_tsv("../processed_data/KD_boolean_matrix_wHead.tsv") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::filter(n_NA < 1)  

#read all Australia Divergent SNP sites, filter sites with >10% missing genotypes
AD_sites <- readr::read_tsv("../processed_data/AD_boolean_matrix_wHead.tsv") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::filter(n_NA < 3) 

#read all Taiwan-Hawaii SNP sites, filter sites with >10% missing genotypes
TH_sites <- readr::read_tsv("../processed_data/TH_boolean_matrix_wHead.tsv") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::filter(n_NA < 10) 

#read all Temperate SNP sites, filter sites with >10% missing genotypes
TEMP_sites <- readr::read_tsv("../processed_data/Temperate_boolean_matrix_wHead.tsv") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::filter(n_NA < 3) 

#read all Taiwan Divergent 1 SNP sites, filter sites with >10% missing genotypes
TD1_sites <- readr::read_tsv("../processed_data/TD1_boolean_matrix_wHead.tsv") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::filter(n_NA < 3) 

#combined shared SNP sites, filter sites with >10% missing genotypes and <1% MAF among Tropical samples
combined_shared <- bind_rows(KDshared  %>% mutate(source = "KD"),
                             ADshared  %>% mutate(source = "AD"),
                             TD1shared %>% mutate(source = "TD1"),
                             TEMPshared %>% mutate(source = "TEMP"),
                             THshared  %>% mutate(source = "TH")) %>%
  dplyr::group_by(across(-source)) %>%  # group by all columns except 'source'
  dplyr::summarise(source = paste(unique(source), collapse = ","),.groups = "drop") %>%
  dplyr::mutate(n_NA = rowSums(is.na(select(., -POS)))) %>%
  dplyr::mutate(n_ones = rowSums(select(., -POS) == 1, na.rm = TRUE)) %>%
  dplyr::mutate(n_zeros = rowSums(select(., -POS) == 0, na.rm = TRUE)) %>%
  dplyr::filter(n_NA < 50) %>% #missing <10%
  dplyr::filter(n_zeros >= 5 & n_ones >= 5) # MAF > 1%


#when deduplicating sites, we want to keep track of what relatedness groups  each site is shared with across different pairwise relatedness group comparisons
#this function collapses relatedness group labels into a comma separated list when de-duplicating sites found in different relatedness group comparisons 
subset_labels <- function(v) {
  v <- sort(unique(v))
  n <- length(v)
  unlist(lapply(seq_len(n), function(k) {
    m <- utils::combn(v, k)
    if (is.null(dim(m))) {
      paste(m, collapse = ",")
    } else {
      apply(m, 2, paste, collapse = ",")
    }
  }))
}

#count the number of shared SNP sites across pairwise and multiple (3-6) relatednessgroup comparisons
sharing <- combined_shared %>%
  dplyr::transmute(POS, 
                   source = paste0("Tropical,", gsub("TEMP", "Temperate", source))) %>%
  dplyr::mutate(tokens = base::strsplit(source, ",", fixed = TRUE)) %>%
  dplyr::mutate(subset = purrr::map(tokens, subset_labels)) %>%
  tidyr::unnest_longer(col = subset, values_to = "subset") %>%
  dplyr::count(subset, name = "num_site") %>%
  dplyr::mutate(groups = stringr::str_count(subset, ",") + 1,
                subset = forcats::as_factor(subset)) %>%
  dplyr::arrange(dplyr::desc(num_site)) %>%
  dplyr::filter(grepl(",",subset) & grepl("Tropical",subset)) %>%
  dplyr::mutate(groups = factor(groups),
                subset = forcats::fct_reorder(subset, num_site, .desc = F))

#visualize the counts of shared SNPs across different relatedness groups comparisons (from 2 to 6 groups compared at a time)
group_joins <- ggplot(data=sharing,aes(x = num_site, y = subset)) +
  geom_col() +
  geom_text(aes(label = num_site),
            hjust = -0.2,              # nudge text right of the bar
            size = 3.2) +
  ggh4x::facet_grid2(
    vars(groups),
    scales = "free_y",   # each facet has its own y (categories) scale
    space  = "free_y"    # facet height scales with number of rows in that facet
  ) +
  labs(x = "Number of shared SNPs",
       y = "") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line()) +
  scale_x_continuous(expand = c(0.01,0))+
  coord_cartesian(xlim = c(0,8.7e4))

#extract filtered SNPs from each relatedness group while keeping track of the shared labels
TH_filt <- combined_shared %>% dplyr::filter(grepl("TH",source))
TEMP_filt <- combined_shared %>% dplyr::filter(grepl("TEMP",source))
TD1_filt <- combined_shared %>% dplyr::filter(grepl("TD1",source))
AD_filt <- combined_shared %>% dplyr::filter(grepl("AD",source))
KD_filt <- combined_shared %>% dplyr::filter(grepl("KD",source))

#Note:
#in the following functions, you will see the term "tlSNPs" invoked often 
#prior to defining C.b. relatedness groups, we referred to these groups as "lineages"
#as such, shared SNPs were initially referred to as trans-lineage SNPs (tlSNPs)
#this collection of functions aims to identify genomic bins that are dense in shared SNPs or show structure in how shared SNPs are distributed across the bin

#for a given genotype matrix of shared SNPs, count the number of shared SNP alternative alleles in every 1 kb genomic bin for each isotype
#statistics about how shared SNP alternative aleles are distributed within the 1 kb bin are also gathered (sparsity (mean inter-SNP distance) and its standard deviation)
count_tlSNPs_per_bin_matrix <- function(bool_df, bins_df) {
  
  #convert to matrix
  pos_col <- grep("^POS", colnames(bool_df))[1]
  pos_vec <- bool_df[[pos_col]]
  pos_split <- str_split_fixed(pos_vec, ":", 2)
  chroms <- pos_split[, 1]
  pos <- as.integer(pos_split[, 2])
  
  mat <- as.matrix(bool_df[, -(pos_col)])
  rownames(mat) <- paste0(chroms, ":", pos)
  colnames(mat) <- colnames(bool_df)[-(pos_col)]
  
  #find ALT==1 entries
  alt_hits <- which(mat == 1, arr.ind = TRUE)
  row_ids <- rownames(mat)[alt_hits[, "row"]]
  col_ids <- colnames(mat)[alt_hits[, "col"]]
  
  #build data.table with CHROM, POS, sample
  split_pos <- data.table::tstrsplit(row_ids, ":", fixed = TRUE)
  dt <- data.table(
    CHROM = split_pos[[1]],
    POS   = as.integer(split_pos[[2]]),
    sample = col_ids
  )
  
  #assign bins and compute count + sparsity score + SD
  dt[, bin_start := floor(POS / 1000) * 1000]
  metrics <- dt[, {
    .N -> count
    sorted_pos <- sort(POS)
    if (count > 1L) {
      gaps <- diff(sorted_pos)
      sparsity <- mean(gaps)
      sparsity_sd <- sd(gaps)
    } else {
      sparsity <- NA_real_
      sparsity_sd <- NA_real_
    }
    .(tlSNP_count = count, sparsity_score = sparsity, sparsity_sd = sparsity_sd)
  }, by = .(CHROM, bin_start, sample)]
  
  #join bin ends
  bins_dt <- as.data.table(bins_df)
  all_samples <- unique(metrics$sample)
  
  #create full bin × sample grid
  bins_expanded <- bins_dt[, .(CHROM, start)][
    , .(sample = all_samples), by = .(CHROM, start)
  ]
  
  bins_expanded <- merge(bins_expanded, bins_dt, by = c("CHROM", "start"), all.x = TRUE)
  
  #merge metrics onto full grid
  result <- merge(
    bins_expanded, metrics,
    by.x = c("CHROM", "start", "sample"),
    by.y = c("CHROM", "bin_start", "sample"),
    all.x = TRUE
  )
  
  #fill NAs
  result[, tlSNP_count := fifelse(is.na(tlSNP_count), 0L, tlSNP_count)]
  # sparsity_score and sparsity_sd remain NA for bins with 0 or 1 SNP
  
  #reshape to wide formats
  wide_counts <- dcast(
    result,
    CHROM + start + end ~ sample,
    value.var = "tlSNP_count",
    fill = 0
  )
  
  wide_sparsity <- dcast(
    result,
    CHROM + start + end ~ sample,
    value.var = "sparsity_score"
  )
  
  wide_sparsity_sd <- dcast(
    result,
    CHROM + start + end ~ sample,
    value.var = "sparsity_sd"
  )
  
  return(list(
    long_result = result,
    wide_counts = wide_counts,
    wide_sparsity = wide_sparsity,
    wide_sparsity_sd = wide_sparsity_sd
  ))
}

#this function filters out genomic bins with shared SNPs that don't meet certain criteria (described below) regarding sparsity and physical distribution
#
filter_tlSNP_data <- function(dtc,scaling_factor) {
  # Work directly on a copy of the input
  dt <- copy(dtc$long_result)
  
  # Calculate combined and z-transformed features
  dt[, combined_sparsity := sparsity_score + sparsity_sd]
  z_mean <- scale(dt$sparsity_score)
  z_sd <- scale(dt$sparsity_sd)
  dt[, z_combined := z_mean + z_sd]
  dt[, z_sd := z_sd]
  dt[, z_mean := z_mean]
  raw <- dt[!is.na(sparsity_score) & !is.na(sparsity_sd)] %>%
    dplyr::filter(tlSNP_count > 3) 
  shift_mean <- -(quantile(raw$z_mean,0.5))
  shift_sd <- -(quantile(raw$z_sd,0.5))
  
  # Apply  filters
  filtered <- raw %>%
    dplyr::filter(ifelse(tlSNP_count < 10 & z_sd < z_mean, FALSE, TRUE)) %>% #1) for bins with less than 10 shared SNPs, remove if the  Z(mean inter-SNP distance) is greater than Z(standard deviation of inter-SNP distance)  (under y=x)
    dplyr::filter(ifelse(z_mean > 0 & z_sd < 0, FALSE, TRUE)) %>% #2) remove bins with Z(mean inter-SNP distance) > 0 and Z(standard deviation of inter-SNP distance)  < 0 - (sparse and unstructured)
    dplyr::filter(ifelse(tlSNP_count < 10 & z_sd+shift_sd < (z_mean+shift_mean) * scaling_factor , FALSE, TRUE)) #3) similar to 1), but we shift from filtering out bins under y=x to bins under y=f(x-s(x)) + s(y) , 
                                                                                                                 #where s(x/y) is the median-shifted Z(mean inter-SNP distance)/Z(standard deviation of inter-SNP distance) 
                                                                                                                 #and f is a scaling factor of the slope (later defined at 2 when the function is invoked).
                                                                                                                 #This third step more aggressively prunes bins with lower shared SNPs counts, by raising the stringency at which a bin with lower
                                                                                                                 #SNP counts is considered to be structured.
  return(list(filtered = filtered, raw = raw))
}

#visualize dotplot of genomic bin Z(SD) and Z(mean) after filtering given the criteria above
plot_filter_diagnostics <- function(dt_filtered, refdt, scaling_factor) {
  #compute percentiles and labels
  percentiles <- quantile(dt_filtered$tlSNP_count, probs = c(0.5, 0.75, 0.95, 0.99, 0.999), na.rm = TRUE)
  percentile_labels <- as.character(percentiles)
  percentile_labels[1] <- paste0("<", percentile_labels[1])
  percentile_labels[length(percentile_labels)] <- paste0(">", percentile_labels[length(percentile_labels)])
  
  #calculate shift values from reference data
  shift_mean <- -quantile(refdt$z_mean, 0.5)
  shift_sd <- -quantile(refdt$z_sd, 0.5)
  
  #define diagonal filter line
  intercept_diag <- scaling_factor * shift_mean - shift_sd
  
  #line definitions
  line_df <- data.frame(
    type = factor(
      c("Z(μ/σ) = 0", "Z(μ/σ) = 0", "Z(μ/σ) = θ_Z(μ/σ)", "Z(μ/σ) = θ_Z(μ/σ)", "Z(σ) = 2(Z(μ) - θ_Z(μ)) + θ_Z(σ)", "Z(σ) = Z(μ)"),
      levels = c("Z(μ/σ) = 0", "Z(μ/σ) = θ_Z(μ/σ)", "Z(σ) = Z(μ)", "Z(σ) = 2(Z(μ) - θ_Z(μ)) + θ_Z(σ)")
    ),
    xintercept = c(0, NA, -shift_mean, NA, NA, NA),
    yintercept = c(NA, 0, NA, -shift_sd, NA, NA),
    slope = c(NA, NA, NA, NA, scaling_factor, 1),
    intercept = c(NA, NA, NA, NA, intercept_diag, 0)
  )
  
  #build plot
  ggplot() +
    geom_point(
      data = dt_filtered,
      aes(x = z_mean, y = z_sd, color = tlSNP_count),
      alpha = 0.6, size = 1
    ) +
    
    geom_vline(
      data = line_df[!is.na(line_df$xintercept) & line_df$type == "Z(μ/σ) = 0", ],
      aes(xintercept = xintercept, linetype = type),
      linewidth = 0.6, color = "black"
    ) +
    geom_vline(
      data = line_df[!is.na(line_df$xintercept) & line_df$type == "Z(μ/σ) = θ_Z(μ/σ)", ],
      aes(xintercept = xintercept, linetype = type),
      linewidth = 0.6, color = "black"
    ) +
    geom_hline(
      data = line_df[!is.na(line_df$yintercept) & line_df$type == "Z(μ/σ) = 0", ],
      aes(yintercept = yintercept, linetype = type),
      linewidth = 0.6, color = "black"
    ) +
    geom_hline(
      data = line_df[!is.na(line_df$yintercept) & line_df$type == "Z(μ/σ) = θ_Z(μ/σ)", ],
      aes(yintercept = yintercept, linetype = type),
      linewidth = 0.6, color = "black"
    ) +
    geom_abline(
      data = line_df[!is.na(line_df$slope) & line_df$type != "Z(σ) = Z(μ)", ],
      aes(slope = slope, intercept = intercept, linetype = type),
      linewidth = 0.6, color = "black"
    ) +
    geom_abline(
      data = line_df[line_df$type == "Z(σ) = Z(μ)", ],
      aes(slope = slope, intercept = intercept, linetype = type),
      linewidth = 0.6, color = "black"
    ) +
    
    scale_linetype_manual(
      name = "",
      values = c(
        "Z(μ/σ) = 0" = "solid",
        "Z(μ/σ) = θ_Z(μ/σ)" = "dashed",
        "Z(σ) = Z(μ)" = "dotted",
        "Z(σ) = 2(Z(μ) - θ_Z(μ)) + θ_Z(σ)" = "dotdash"
      )
    ) +
    
    scale_color_viridis_c(
      option = "plasma",
      limits = range(percentiles),
      breaks = percentiles,
      labels = percentile_labels,
      oob = scales::squish,
      name = "mean\ntlSNP count",
      guide = guide_colorbar(
        ticks = TRUE,
        ticks.colour = "white",
        ticks.linewidth = 0.4,
        frame.colour = "black",
        label.hjust = 1
      )
    ) +
    labs(x = "", y = "") +
    theme_bw()+ 
    theme(legend.justification.right = "top")
}



#visualize dotplot of genomic bin Z(SD) and Z(mean) before filtering 
plot_filter_diagnostics_raw <- function(refdt) {
  #compute percentiles and labels
  percentiles <- quantile(refdt$tlSNP_count, probs = c(0.5, 0.75, 0.95, 0.99, 0.999), na.rm = TRUE)
  percentile_labels <- as.character(percentiles)
  percentile_labels[1] <- paste0("<", percentile_labels[1])
  percentile_labels[length(percentile_labels)] <- paste0(">", percentile_labels[length(percentile_labels)])
  
  ggplot() +
    geom_point(data=refdt, aes(x = z_mean, y = z_sd, color = tlSNP_count), alpha = 0.6, size = 1) +
    scale_color_viridis_c(
      option = "plasma",
      limits = range(percentiles),
      breaks = percentiles,
      labels = percentile_labels,
      oob = scales::squish,
      name = "mean\ntlSNP count",
      guide = guide_colorbar(
        ticks = TRUE,
        ticks.colour = "white",
        ticks.linewidth = 0.4,
        frame.colour = "black",
        label.hjust = 1
      )
    ) +
    labs(x = "", y = "") +
    theme_bw() + 
    theme(legend.justification.right = c(0,1))
}


#visualize heatmap of the mean number of shared SNPs across different values of Z(SD) and Z(mean)
plot_z_heatmap <- function(dt, refdt, bin_width = 0.1, title = "") {
  dt_copy <- copy(dt)
  refdt_copy <- copy(refdt)
  
  #define bin edges from refdt
  min_zsd <- floor(min(refdt_copy$z_sd, na.rm = TRUE) / bin_width) * bin_width
  max_zsd <- ceiling(max(refdt_copy$z_sd, na.rm = TRUE) / bin_width) * bin_width
  breaks_zsd <- seq(min_zsd, max_zsd, by = bin_width)
  
  min_zmean <- floor(min(refdt_copy$z_mean, na.rm = TRUE) / bin_width) * bin_width
  max_zmean <- ceiling(max(refdt_copy$z_mean, na.rm = TRUE) / bin_width) * bin_width
  breaks_zmean <- seq(min_zmean, max_zmean, by = bin_width)
  
  #bin midpoints for axis limits
  bin_mid <- function(breaks) head(breaks, -1) + diff(breaks)/2
  zsd_midpoints <- bin_mid(breaks_zsd)
  zmean_midpoints <- bin_mid(breaks_zmean)
  
  #bin the data
  dt_copy[, zsd_bin := cut(z_sd, breaks = breaks_zsd, include.lowest = TRUE)]
  dt_copy[, zmean_bin := cut(z_mean, breaks = breaks_zmean, include.lowest = TRUE)]
  
  #map bins to midpoints
  zsd_levels <- levels(dt_copy$zsd_bin)
  zmean_levels <- levels(dt_copy$zmean_bin)
  zsd_mid_map <- setNames(zsd_midpoints, zsd_levels)
  zmean_mid_map <- setNames(zmean_midpoints, zmean_levels)
  
  dt_copy[, zsd_mid := zsd_mid_map[as.character(zsd_bin)]]
  dt_copy[, zmean_mid := zmean_mid_map[as.character(zmean_bin)]]
  
  #count observations per bin (only observed ones)
  heatmap_data <- dt_copy[, .N, by = .(zmean_mid, zsd_mid)]
  
  N_percentiles <- quantile(heatmap_data$N, probs = c(0.5, 0.75, 0.95,0.99,0.999), na.rm = TRUE)
  percentile_labels <- as.character(round(N_percentiles))
  percentile_labels[1] <- paste0("<", percentile_labels[1])
  percentile_labels[length(percentile_labels)] <- paste0(">", percentile_labels[length(percentile_labels)])
  
  #plot with fixed x/y limits based on reference breaks
  ggplot(heatmap_data, aes(x = zmean_mid, y = zsd_mid, fill = N)) +
    geom_tile(color = "white") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_x_continuous(
      limits = range(zmean_midpoints),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = range(zsd_midpoints),
      expand = c(0, 0)
    ) +
    scale_fill_viridis_c(
      option = "plasma",
      limits = range(N_percentiles),
      breaks = N_percentiles,
      labels = percentile_labels,
      oob = scales::squish,
      name = "number\nof bins",
      guide = guide_colorbar(
        ticks = TRUE,
        ticks.colour = "white",
        ticks.linewidth = 0.4,
        frame.colour = "black",
        label.hjust = 1
      )
    ) +
    labs(
      title = title,
      x = "",
      y = "",
      fill = "bin count"
    ) +
    theme_bw()
}

#visualize heatmap of the number of genomic bins at different values of Z(SD) and Z(mean)
plot_z_heatmap_count <- function(dt, refdt, bin_width = 0.1, title = "") {
  dt_copy <- copy(dt)
  refdt_copy <- copy(refdt)
  
  #define bin edges from refdt
  min_zsd <- floor(min(refdt_copy$z_sd, na.rm = TRUE) / bin_width) * bin_width
  max_zsd <- ceiling(max(refdt_copy$z_sd, na.rm = TRUE) / bin_width) * bin_width
  breaks_zsd <- seq(min_zsd, max_zsd, by = bin_width)
  
  min_zmean <- floor(min(refdt_copy$z_mean, na.rm = TRUE) / bin_width) * bin_width
  max_zmean <- ceiling(max(refdt_copy$z_mean, na.rm = TRUE) / bin_width) * bin_width
  breaks_zmean <- seq(min_zmean, max_zmean, by = bin_width)
  
  #bin midpoints
  bin_mid <- function(breaks) head(breaks, -1) + diff(breaks)/2
  zsd_midpoints <- bin_mid(breaks_zsd)
  zmean_midpoints <- bin_mid(breaks_zmean)
  
  #bin the data
  dt_copy[, zsd_bin := cut(z_sd, breaks = breaks_zsd, include.lowest = TRUE)]
  dt_copy[, zmean_bin := cut(z_mean, breaks = breaks_zmean, include.lowest = TRUE)]
  
  #map bin factors to midpoints
  zsd_levels <- levels(dt_copy$zsd_bin)
  zmean_levels <- levels(dt_copy$zmean_bin)
  zsd_mid_map <- setNames(zsd_midpoints, zsd_levels)
  zmean_mid_map <- setNames(zmean_midpoints, zmean_levels)
  
  dt_copy[, zsd_mid := zsd_mid_map[as.character(zsd_bin)]]
  dt_copy[, zmean_mid := zmean_mid_map[as.character(zmean_bin)]]
  
  #compute mean tlSNP per bin
  heatmap_data <- dt_copy[, .(mean_tlSNP = mean(tlSNP_count, na.rm = TRUE)), by = .(zmean_mid, zsd_mid)]
  
  #compute reference percentiles from ALL tlSNP values
  tlSNP_percentiles <- quantile(heatmap_data$mean_tlSNP, probs = c(0.75, 0.95,0.99,0.999), na.rm = TRUE)
  percentile_labels <- as.character(round(tlSNP_percentiles))
  percentile_labels[1] <- paste0("<", percentile_labels[1])
  percentile_labels[length(percentile_labels)] <- paste0(">", percentile_labels[length(percentile_labels)])
  #scale_limits <- range(tlSNP_percentiles)
  
  #plot
  ggplot(heatmap_data, aes(x = zmean_mid, y = zsd_mid, fill = mean_tlSNP)) +
    geom_tile(color = "white") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_x_continuous(limits = range(zmean_midpoints), expand = c(0, 0)) +
    scale_y_continuous(limits = range(zsd_midpoints), expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "plasma",
      limits = range(tlSNP_percentiles),
      breaks = tlSNP_percentiles,
      labels = percentile_labels,
      oob = scales::squish,
      name = "mean shared\nSNP count",
      guide = guide_colorbar(
        ticks = TRUE,
        ticks.colour = "white",
        ticks.linewidth = 0.4,
        frame.colour = "black",
        label.hjust = 1
      )
    ) +
    labs(
      title = title,
      x = "",
      y = "",
      fill = "bin count"
    ) +
    theme_bw()
}

#combine plots and generate plot list
plot_all_diagnostics <- function(filtered_data, reference_data, scaling_factor = 2) {
  # Axis labels
  x_lab <- ggdraw() + draw_label("Z(mean inter-SNP distance)")
  y_lab <- ggdraw() + draw_label("Z(SD of inter-SNP distance)", angle = 90)
  
  # Plot diagnostics
  dig_filtered <- plot_filter_diagnostics(filtered_data, reference_data, scaling_factor)
  dig_raw <- plot_filter_diagnostics_raw(reference_data)
  
  combined_dig <- plot_grid(
    dig_raw, dig_filtered, NULL,
    ncol = 3,
    align = "hv",
    labels = c("", "b"),
    rel_widths = c(1, 1, 0.02)
  )
  
  dig_plot <- plot_grid(
    plot_grid(NULL, y_lab, NULL, ncol = 1, rel_heights = c(0.01, 1, 0.01)),
    plot_grid(combined_dig, x_lab, NULL, ncol = 1, rel_heights = c(1, 0.01, 0.01)),
    ncol = 2,
    rel_widths = c(0.02, 1),
    labels = c("a", "")
  )
  
  #heatmaps for bin counts
  hm_ref <- plot_z_heatmap(reference_data, reference_data)
  hm_filtered <- plot_z_heatmap(filtered_data, reference_data)
  
  combined_binct <- plot_grid(hm_ref, hm_filtered, ncol = 2, align = "hv", labels = c("", "b"))
  
  binct_plot <- plot_grid(
    plot_grid(NULL, y_lab, NULL, ncol = 1, rel_heights = c(0.01, 1, 0.01)),
    plot_grid(combined_binct, x_lab, NULL, ncol = 1, rel_heights = c(1, 0.01, 0.01)),
    ncol = 2,
    rel_widths = c(0.02, 1),
    labels = c("a", "")
  )
  
  #heatmaps for mean tlSNP per bin
  hm_count_ref <- plot_z_heatmap_count(reference_data, reference_data)
  hm_count_filtered <- plot_z_heatmap_count(filtered_data, reference_data)
  
  combined_mean <- plot_grid(hm_count_ref, hm_count_filtered, ncol = 2, align = "hv", labels = c("", "b"))
  
  bin_meantl <- plot_grid(
    plot_grid(NULL, y_lab, NULL, ncol = 1, rel_heights = c(0.01, 1, 0.01)),
    plot_grid(combined_mean, x_lab, NULL, ncol = 1, rel_heights = c(1, 0.01, 0.01)),
    ncol = 2,
    rel_widths = c(0.02, 1),
    labels = c("a", "")
  )
  
  return(list(
    combined_dig = dig_plot,
    bin_meantl = bin_meantl,
    binct_plot = binct_plot
  ))
}

#invoke the functions above for each pairiwse relatedness group comparison
#all comparisons are relative to Tropical relatedness group
scaling_factor=2
KD_counts <- count_tlSNPs_per_bin_matrix(KD_filt %>% dplyr::select(-source,-n_NA,-n_ones,-n_zeros), bins)
KD_result <- filter_tlSNP_data(KD_counts,scaling_factor)
KD_filtered <- KD_result$filtered
KD_refpt <- KD_result$raw
KD_plots <-plot_all_diagnostics(KD_filtered,KD_refpt,scaling_factor)

AD_counts <- count_tlSNPs_per_bin_matrix(AD_filt %>% dplyr::select(-source,-n_NA,-n_ones,-n_zeros), bins)
AD_result <- filter_tlSNP_data(AD_counts,scaling_factor)
AD_filtered <- AD_result$filtered
AD_refpt <- AD_result$raw
AD_plots <-plot_all_diagnostics(AD_filtered,AD_refpt,scaling_factor)

TD1_counts   <- count_tlSNPs_per_bin_matrix(TD1_filt %>% dplyr::select(-source,-n_NA,-n_ones,-n_zeros), bins)
TD1_result <- filter_tlSNP_data(TD1_counts,scaling_factor)
TD1_filtered <- TD1_result$filtered
TD1_refpt <- TD1_result$raw
TD1_plots <-plot_all_diagnostics(TD1_filtered,TD1_refpt,scaling_factor)

TMP_counts   <- count_tlSNPs_per_bin_matrix(TEMP_filt %>% dplyr::select(-source,-n_NA,-n_ones,-n_zeros), bins)
TMP_result <- filter_tlSNP_data(TMP_counts,scaling_factor)
TMP_filtered <- TMP_result$filtered
TMP_refpt <- TMP_result$raw
TMP_plots <-plot_all_diagnostics(TMP_filtered,TMP_refpt,scaling_factor)

TH_counts   <- count_tlSNPs_per_bin_matrix(TH_filt %>% dplyr::select(-source,-n_NA,-n_ones,-n_zeros), bins)
TH_result <- filter_tlSNP_data(TH_counts,scaling_factor)
TH_filtered <- TH_result$filtered
TH_refpt <- TH_result$raw
TH_plots <-plot_all_diagnostics(TH_filtered,TH_refpt,scaling_factor)

############################################# shared SNP and HDR overlap ###########################################

#combined shared contains isotype-specific genotypes for shared SNP sites, and keeps track of which relatedness group each site is shared among
#split POS, reshape to long
dt <- as.data.table(combined_shared %>% dplyr::select(-n_NA,-n_ones,-n_zeros))
dt[, c("region_chrom","site") := data.table::tstrsplit(POS, ":", fixed = TRUE)]
dt[, site := as.integer(site)]

#melt to long, STRAIN = sample column name, genotype = value
long <- melt(
  dt,
  id.vars = c("POS","region_chrom","site","source"),
  variable.name = "STRAIN",
  value.name = "genotype"
)

site_ranges <- long[, .(STRAIN, POS, genotype,source,
                        region_chrom, start = site, end = site)] #isotype-specific genotypes

unique_sites <- site_ranges %>% dplyr::select(-STRAIN,-genotype) %>% dplyr::distinct(POS,.keep_all = T) #remove STRAIN and genotype, just keep unique sites


hdrs_dt <- as.data.table(hdrs)[, .(STRAIN, region_chrom = CHROM,
                                   region_start = minStart,
                                   region_end = maxEnd)] #hdrs

bins_dt <-as.data.table(bins)[, .(region_chrom = CHROM,
                                  bin_start = start,
                                  bin_end = end)] #genomic bins

#set keys for foverlaps()
setkey(site_ranges, STRAIN, region_chrom, start, end)
setkey(unique_sites, region_chrom, start, end)
setkey(hdrs_dt,   STRAIN, region_chrom, region_start, region_end) 
setkey(bins_dt, region_chrom, bin_start, bin_end)

setnames(hdrs_dt, c("region_start","region_end"), c("start","end"))
setnames(bins_dt, c("bin_start","bin_end"), c("start","end"))
ov <- foverlaps(site_ranges, hdrs_dt, type = "any") #find overlaps between shared SNPs and HDRs in each strain, to count the number of relatedness groups that share SNPs in each hyper-divergent region
ov_binned <- foverlaps(unique_sites, bins_dt, type = "any") #find genomic bins that overlap with shared SNPs, to estimate the proportion of the reference genome that contains shared SNPs for each pairwise comparison

ov_long <- ov_binned[, .(source = trimws(unlist(strsplit(source, ",")))),
                  by = .(region_chrom, start, end)] #generate an individual row for each relatedness group which whom the SNP is shared with

ov_long_unique <- unique(ov_long) #some bins contain several different SNPs shared with the same relatedness group
                                  #we dedup, since we only want to count that bin once to have an accurate span estimate when summing the count of bins

src_counts <- ov_long_unique[, .N, by = source][order(-N)] #sum the number of bins for each pairwise comparison

#plot the proportion of the genome with shared SNPs
pct_shared <- ggplot(src_counts %>% dplyr::mutate(source=ifelse(source=="TEMP","Temperate",source)),
       aes(x = forcats::fct_rev(forcats::fct_reorder(source, N, .desc = TRUE)),
           y = (N / nrow(bins)) * 100)) +
  geom_col(width = 0.7) +
  geom_text(aes(
    label = paste0(scales::comma(round((N / nrow(bins) * 100), 1)), "%") # nrow(bins) is the number of 1kb genomic bins, which is effectively the genome size in kb
  ),
  hjust = -0.2, size = 3.2) +
  labs(x = "",
       y = "Percent genome span with SNPs\nshared with Tropical group") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    plot.margin = margin(10,5,5,5)) +
  scale_y_continuous(expand = c(0.01, 0.1), limits = c(0, 23)) +
  coord_flip()

#count the total number of SNPs in each relatedness group 
df_counts <- tibble::tibble(
  group = c("Tropical","KD","AD","TH","Temperate","TD1"),
  n     = c(nrow(Trop_sites),
            nrow(KD_sites),
            nrow(AD_sites),
            nrow(TH_sites),
            nrow(TEMP_sites),
            nrow(TD1_sites))
) %>%
  dplyr::mutate(group = forcats::fct_reorder(group, n, .desc = TRUE))

#side bar chart of SNP counts
side_counts <- ggplot(df_counts, aes(x = n, y = forcats::fct_rev(group))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15, size = 3.2) +
  labs(x = "Number of SNPs", y = "Relatedness\ngroup") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line()
  ) +
  scale_x_continuous(expand = c(0.01, 0.15),limits = c(0,nrow(Trop_sites)+2.5e5))

#combined plot of SNP counts, shared SNP counts, and proportion of the genome with shared SNPs across different relatedness group comparisons (Fig S23)
bc_pan <-cowplot::plot_grid(side_counts,pct_shared,nrow=1,align = "h",axis="tb",labels=c("a","b"))
abc_pan <- cowplot::plot_grid(bc_pan,group_joins,nrow=2,rel_heights = c(1,2.5),labels=c("","c"))
ggsave(plot = abc_pan, filename = "../figures/FigureS23_shared_snp.png",width = 7,height = 6,dpi = 600,device = 'png',bg = "white")


#shared SNPs in HDRs within and outside of coding sequences
#we want to identify how many different relatedness groups SNPs are shared with in HDRs
ov_dt  <- as.data.table(ov_df)
cds_dt <- as.data.table(cds)

ov_dt[, idx := .I]
ov_dt[, site := as.integer(site)]
ov_dt[, `:=`(start = site, end = site)]

cds_dt[, `:=`(start = as.integer(start), end = as.integer(end))]

setkey(ov_dt,  CHROM, start, end)
setkey(cds_dt, chrom, start, end)

cds_ov_dt <- foverlaps(
  x    = ov_dt,
  y    = cds_dt,
  by.x = c("CHROM","start","end"),
  by.y = c("chrom","start","end")
) # overlaps with coding sequences

cds_hits <- as.data.frame(cds_ov_dt) %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::group_by(STRAIN,POS) %>%
  dplyr::distinct(POS,.keep_all = T) %>%
  dplyr::ungroup() #hits with CDS

cds_props <- cds_hits %>% 
  dplyr::filter(genotype==1) %>%  #alt alleles only
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(total_alt_sites=sum(!is.na(genotype)), overlapping_alt_sites=sum(!is.na(region_start))) %>%
  dplyr::mutate(prop_overlap=overlapping_alt_sites/total_alt_sites) 

#get number of relatedness group CDS SNPs alt alleles contained in each HDR
cds_combos <- cds_hits %>%
  dplyr::filter(!is.na(region_start)) %>%
  dplyr::group_by(STRAIN, CHROM, region_start, region_end) %>%
  dplyr::summarise(
    # split all sources in the group on ",", flatten, trim, and count unique tokens
    num_lin = dplyr::n_distinct(
      trimws(
        unlist(base::strsplit(stats::na.omit(source), ",", fixed = TRUE))
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    midpoint = region_start + ((region_end - region_start) / 2),
    size = region_end - region_start
  )

#get number of relatedness group SNPs alt alleles contained in each HDR
all_combos <- ov_df %>%
  dplyr::filter(!is.na(region_start))  %>%
  dplyr::group_by(STRAIN, CHROM, region_start, region_end) %>%
  dplyr::summarise(
    # split all sources in the group on ",", flatten, trim, and count unique tokens
    num_lin = dplyr::n_distinct(
      trimws(
        unlist(base::strsplit(stats::na.omit(source), ",", fixed = TRUE))
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    midpoint = region_start + ((region_end - region_start) / 2),
    size = region_end - region_start
  )

#plot
plot_combos <- rbind(all_combos %>% dplyr::mutate(class="All shared SNPs"),
                     cds_combos %>% dplyr::mutate(class="Shared SNPs in\nprotein-coding sequences"))
combos_plt <- ggplot(plot_combos) + 
  geom_jitter(aes(y=size/1e3,x=num_lin,group=num_lin),alpha=0.4,shape=21,color="forestgreen",height = 1) + 
  geom_boxplot(aes(y=size/1e3,x=num_lin,group=num_lin),outliers = F,color="black",fill=NA) + theme_bw() + 
  xlab("Number of relatedness groups") + 
  ylab("Size of hyper-divergent region (kb)") +
  facet_wrap(~class,ncol=2) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())+
  scale_y_continuous(expand = c(0.01,0),limits = c(0,220))

ggsave(combos_plt,filename = "../figures/FigureS24_RGcount_perRegion.png",device = "png",dpi = 600,units = "in",width = 7.5,height = 7.5,bg="white")

lineages <- readr::read_tsv("../processed_data/isotype_byLineage_GeoLocAdmCol_20250909.tsv") %>%
  dplyr::mutate(sublineage_color=ifelse(Sublineage=="TC","#ff0000",sublineage_color)) %>%
  dplyr::mutate(Sublineage=ifelse(Sublineage=="TC","TT",Sublineage)) 

geno_only <- ADshared %>% select(-POS)
na_props <- colMeans(is.na(geno_only))
keep_cols <- names(na_props[na_props <= 0.10])
filtered_df_clean <- ADshared %>% select(POS, all_of(keep_cols))

filtered_df2 <- filtered_df_clean %>%
  tibble::column_to_rownames(var = "POS")

chrom_order <- c("I", "II", "III", "IV", "V", "X")

# Extract CHROM and POS from row names
filtered_df2_sorted <- filtered_df2 %>%
  tibble::rownames_to_column(var = "variant") %>%
  separate(variant, into = c("CHROM", "POS"), sep = ":", remove = FALSE) %>%
  mutate(
    CHROM = factor(CHROM, levels = chrom_order),
    POS = as.numeric(POS)
  ) %>%
  arrange(CHROM, POS) %>%
  mutate(variant = paste(CHROM, POS, sep = ":")) %>%
  select(-CHROM, -POS) %>%
  tibble::column_to_rownames(var = "variant") 

filtered_df2_t <- as.data.frame(t(filtered_df2_sorted))

heat_mat_num <- apply(filtered_df2_t, 2, function(x) {
  x <- as.character(x)
  x[x %in% c("NA", "", ".")] <- NA
  as.numeric(x)
})
heat_mat_num <- as.matrix(heat_mat_num)
rownames(heat_mat_num) <- rownames(filtered_df2_t)


strain_levels <- rownames(heat_mat_num)[rownames(heat_mat_num) %in% unique((hdrs)$STRAIN)] #samples_plot[grepl("BRC",samples_plot)]

sublineage_order <- c("TT", "TI", "TA")
geo_order <- function(x) {
  # Push "Globalized" and "unknown" to the end
  x <- ifelse(x %in% c("Globalized", "unknown"), paste0("zzz_", x), x)
  return(x)
}

ordered_strains <- lineages %>%
  dplyr::filter(isotype %in% strain_levels,
         Lineage=="Tropical") %>%
  dplyr::mutate(
    Sublineage = factor(Sublineage, levels = sublineage_order),
    geo_mod = geo_order(geo)
  ) %>%
  dplyr::arrange(Sublineage, geo_mod) %>%
  dplyr::pull(isotype)


# Subset and reorder matrix by strain_levels
mat_sub <- heat_mat_num[ordered_strains, , drop = FALSE]
hit_all <- which(mat_sub == 1 | is.na(mat_sub), arr.ind = TRUE)

# Construct long format data frame with genotypes
df_points_all <- tibble(
  sample   = rownames(mat_sub)[hit_all[, "row"]],
  var      = colnames(mat_sub)[hit_all[, "col"]],
  genotype = mat_sub[hit_all]  # Either 1 or NA
) %>%
  separate_wider_delim(var, ":", names = c("chrom", "pos")) %>%
  dplyr::mutate(
    pos    = as.numeric(pos),
    x      = pos / 1e6,
    sample = factor(sample, levels = ordered_strains)
  ) %>%
  dplyr::mutate(genotype=ifelse(is.na(genotype),"MISS","ALT"))


bin_ranges <- bins %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(xmin = min(start), xmax = max(end), .groups = "drop")

get_xlims_mb <- function(chr) {
  r <- bin_ranges %>% filter(CHROM == chr)
  if (!nrow(r)) stop("No bin range for chromosome ", chr, " in `bins`.")
  c(r$xmin, r$xmax) / 1e6
}


make_joint_track <- function(chr,
                             enable_var=F,
                             kbsp=0.01,
                             angle=0,
                             hjust=0.5,
                             hdr_band = 0.45,
                             hdr_alpha = 0.35,
                             point_size = 0.2,
                             point_alpha = 1,
                             show_y_labels = F) {
  
  hdrs_chr <- hdrs %>%
    filter(CHROM == chr, STRAIN %in% ordered_strains) %>%
    mutate(
      y    = match(STRAIN, ordered_strains),
      xmin = minStart / 1e6,
      xmax = maxEnd   / 1e6,
      ymin = y - hdr_band,
      ymax = y + hdr_band
    )
  
  KD_tlsnpbin <- KD_filtered %>% dplyr::filter(CHROM==chr, sample %in% ordered_strains) %>%
    dplyr::mutate(y = match(sample, ordered_strains),
                  xmin = start / 1e6,
                  xmax = end   / 1e6,
                  ymin = y - hdr_band,
                  ymax = y + hdr_band)
  
  AD_tlsnpbin <- AD_filtered %>% dplyr::filter(CHROM==chr, sample %in% ordered_strains) %>%
    dplyr::mutate(y = match(sample, ordered_strains),
                  xmin = start / 1e6,
                  xmax = end   / 1e6,
                  ymin = y - hdr_band,
                  ymax = y + hdr_band)
  
  TD1_tlsnpbin <- TD1_filtered %>% dplyr::filter(CHROM==chr, sample %in% ordered_strains) %>%
    dplyr::mutate(y = match(sample, ordered_strains),
                  xmin = start / 1e6,
                  xmax = end   / 1e6,
                  ymin = y - hdr_band,
                  ymax = y + hdr_band)
  
  TMP_tlsnpbin <- TMP_filtered %>% dplyr::filter(CHROM==chr, sample %in% ordered_strains) %>%
    dplyr::mutate(y = match(sample, ordered_strains),
                  xmin = start / 1e6,
                  xmax = end   / 1e6,
                  ymin = y - hdr_band,
                  ymax = y + hdr_band)
  
  TH_tlsnpbin <- TH_filtered %>% dplyr::filter(CHROM==chr, sample %in% ordered_strains) %>%
    dplyr::mutate(y = match(sample, ordered_strains),
                  xmin = start / 1e6,
                  xmax = end   / 1e6,
                  ymin = y - hdr_band,
                  ymax = y + hdr_band)
  
  df_points_chr <- df_points_all %>%
    filter(chrom == chr,genotype=="ALT") %>%
    mutate(y = match(sample, ordered_strains))
  
  
  xlims_mb <- get_xlims_mb(chr)
  half_breaks <- seq(
    from = floor(xlims_mb[1] / 0.5) * 0.5,
    to   = ceiling(xlims_mb[2] / 0.5) * 0.5,
    by   = kbsp
  )
  
  ylabs <- if (show_y_labels) ordered_strains else rep("", length(ordered_strains))
  
  if (enable_var==T) {
    ggplot(hdrs_chr) +
      geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "HDR"),
        color = NA, alpha = 1
      ) +
      geom_point(data=df_points_chr, aes(x = x, y = y), color="black", size = point_size, alpha = point_alpha, na.rm = TRUE) +
      geom_rect(
        data=AD_tlsnpbin,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "AD clusters"),
        color = NA, alpha = 0.5
      ) +
    scale_x_continuous(limits = xlims_mb,  breaks = half_breaks, expand = c(0,0)) +
      scale_y_continuous(
        breaks = seq_along(ordered_strains),
        labels = ylabs,
        expand = c(0,0)
      ) +
      scale_fill_manual(
        values = c(
          "HDR"         = "grey50",   
          "AD clusters" = "#FF0000"  
          #"KD clusters" = "#014421"   
          # "TD1 clusters" = "#E6A832",  
          # "Temperate clusters" = "#3C5488",
          # "TH clusters" = "#c700ff"
        ),
        breaks = c("HDR","AD clusters")
      ) +
      labs(x = paste0("Genomic position on chr ", chr, " (Mb)"),
           y = NULL,
           fill="Genomic\nsegments") +        # suppress x label (bottom track shows it)
      theme_bw() +
      theme(
        axis.text.y  = if (show_y_labels) element_text(size = 6) else element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "right",
        legend.justification = "left",
        legend.title.align = 0,
        legend.text.align = 0
      ) + 
      guides(fill = guide_legend(override.aes = list(alpha = 0.6)))
    
  } else {
    
    
    ggplot(hdrs_chr) +
      geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "HDR"),
        color = NA, alpha = 1
      ) +
      geom_rect(
        data=TMP_tlsnpbin,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "TH clusters"),
        color = NA, alpha = 0.5
      ) +
      geom_rect(
        data=KD_tlsnpbin,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "KD clusters"),
        color = NA, alpha = 0.5
      ) +
      geom_rect(
        data=AD_tlsnpbin,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "AD clusters"),
        color = NA, alpha = 0.5
      ) +
      geom_rect(
        data=TD1_tlsnpbin,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "TD1 clusters"),
        color = NA, alpha = 0.5
      ) +
      geom_rect(
        data=TMP_tlsnpbin,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = "Temperate clusters"),
        color = NA, alpha = 0.5
      ) +
      
      
      #
      scale_x_continuous(limits = xlims_mb,  breaks = half_breaks, expand = c(0,0)) +
      scale_y_continuous(
        breaks = seq_along(ordered_strains),
        labels = ylabs,
        expand = c(0,0)
      ) +
      scale_fill_manual(
        values = c(
          "HDR"         = "grey50",   
          "AD clusters" = "#FF0000",  
          "KD clusters" = "#014421",   
          "TD1 clusters" = "#E6A832",  
          "Temperate clusters" = "#3C5488",
          "TH clusters" = "#c700ff"
        ),
        breaks = c("HDR","AD clusters","KD clusters","TD1 clusters", "TH clusters","Temperate clusters")
      ) +
      labs(x = paste0("Genomic position on chr ", chr, " (Mb)"),
           y = NULL,
           fill="Genomic\nsegments") +        # suppress x label (bottom track shows it)
      theme_bw() +
      theme(
        axis.text.y  = if (show_y_labels) element_text(size = 6) else element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "right",
        legend.justification = "left",
        legend.title.align = 0,
        legend.text.align = 0,
        axis.text.x = element_text(angle=angle,hjust=hjust)
      ) + 
      guides(fill = guide_legend(override.aes = list(alpha = 0.6)))  #+
    #scale_color_manual(values=c("MISS"="black","ALT"="blue"))
  }
}



## 5. Stitch HDR + genotype tracks 
plot_chr_stitched <- function(chr, st, en, var, kbsp,
                              angle=0,
                              hjust=0.5,# relative height HDR : points
                              show_y_labels=F) {
  
  ylabs <- if (show_y_labels) ordered_strains else rep("", length(ordered_strains))
  
  lineage_bar <- lineages %>%
    dplyr::filter(isotype %in% ordered_strains) %>%
    dplyr::mutate(
      y = base::match(isotype, ordered_strains),
      ymin = y - 0.5,
      ymax = y + 0.5
    )
  
  p_geo_bar <- ggplot() +
    geom_rect(data = lineage_bar,
              aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = geo),
              color = NA) +
    scale_fill_manual(
      name = "Geography",
      values = setNames(unique(lineage_bar$geo_color), unique(lineage_bar$geo))
    ) +
    scale_y_continuous(
      breaks = seq_along(ordered_strains),
      labels = ylabs,
      expand = c(0,0)
    ) + 
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      legend.position = "none",  # hide legend here
      axis.text.y=element_blank(),
      plot.margin = margin(10, 1, 0, 0)
    )
  
  # Extract geo legend
  p_geo_legend <- suppressWarnings(get_legend(
    p_geo_bar +
      theme(
        legend.position = "right",
        legend.justification = "left",
        legend.title.align = 0,
        legend.text.align = 0
      )))
  
  p_join <- make_joint_track(chr, var, kbsp, angle,hjust, show_y_labels = show_y_labels)
  lims_x = c(st,en)
  p_bool_legend <- suppressWarnings(get_legend(p_join + theme(legend.margin = margin(10,0,0,0))))
  
  
  legend_col <- cowplot::plot_grid(cowplot::ggdraw(p_bool_legend),
                                   #cowplot::ggdraw(p_sublineage_legend),
                                   cowplot::ggdraw(p_geo_legend),
                                   NULL, ncol = 1, align = "v",axis = "lr", rel_heights = c(0.7,1.1,0.25))  
  
  
  pmain <- cowplot::plot_grid(#p_sublineage_bar,
                              p_geo_bar,
                              p_join + theme(legend.position = "none", plot.margin=margin(10,10,0,0))+ coord_cartesian(xlim = lims_x),
                              
                              rel_widths = c(0.5,10),ncol=3,align = "h",axis = "rl")
  
  cowplot::plot_grid(NULL,pmain,legend_col,rel_widths = c(0.1,4,1.2),ncol=3)
}


plt1 <- plot_chr_stitched("IV",4.37,4.41,T,0.01)  
ggsave(plt1,filename = "../figures/FigureS25_IV_4.37-4.41_tlSNP_blockAD_wvar_20251021.png",device = "png",dpi = 600,units = "in",width = 7.5,height = 7.5,bg="white")

plt2 <- plot_chr_stitched("V",4.275,4.377,F,0.01) 
ggsave(plt2,filename = "../figures/FigureS26_V_4.27-4.38_tlSNP_blockAD_20251021.png",device = "png",dpi = 600,units = "in",width = 7.5,height = 7.5,bg="white")

plt3 <- plot_chr_stitched("II",12.92,13.02,F,0.01) 
ggsave(plt3,filename = "../figures/FigureS27_II_12.92-13.02_tlSNP_blockAD_20251021.png",device = "png",dpi = 600,units = "in",width = 7.5,height = 7.5,bg="white")

plt4 <- plot_chr_stitched("V",15.91,15.95,F,0.01)
ggsave(plt4,filename = "../figures/FigureS28_V_15.91-15.95_tlSNP_block_20251021.png",device = "png",dpi = 600,units = "in",width = 7.5,height = 7.5,bg="white")

plt5 <- plot_chr_stitched("I",11.88,12.1,F,0.01,45,1) 
ggsave(plt5,filename = "../figures/FigureS29_I_11.88-12.10_tlSNP_block_20251021.png",device = "png",dpi = 600,units = "in",width = 7.5,height = 7.5,bg="white")
