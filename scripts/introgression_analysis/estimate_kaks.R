library(seqinr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

hdr_groups <- readLines("../../processed_data/tree_topology/aa_divergence/hdr_ogs_all_tropical.txt")
non_hdr_groups <- readLines("../../processed_data/tree_topology/aa_divergence/non_hdr_ogs_all_tropical.txt")
aln_dir <- "../../processed_data/tree_topology/aa_divergence/codon_aware_msa_wgaps"
no_cre_files <-  readLines("../../processed_data/tree_topology/aa_divergence/msa/no_cre_files.txt")
no_cre_files <- gsub(no_cre_files,pattern = ".no_cre.aligned.fa",replacement = "")

hdrs <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_allStrain_5kbclust_20260506.tsv")

gff <- readr::read_tsv("../../processed_data/gene_diversity/c_briggsae.QX1410_20250929.csq.gff", col_names = c("seqid","source","type","start","end","score","strand","phase","attributes")) %>%
  dplyr::filter(seqid!="MtDNA")
qxtran <- gff %>% dplyr::filter(type=="mRNA")

ogs <- readr::read_tsv("../../processed_data/tree_topology/prot_TropicalCB_CN_CR_wDRG_39taxa/OrthoFinder/Results_CBCNCR_39_famsa/Orthogroups/Orthogroups.tsv")
ogs_colnames <- colnames(ogs)
colnames(ogs) <- sub("\\.longest\\.prot\\.otu$", "", ogs_colnames)
nigoni <- c("EG5268","JU2484","JU2617","YR106","ZF1220","ECA2852","ECA2857","JU1418","JU1419","JU1420","JU1422","NIC2143","NIC2150","NIC2152","VX151","VX153") 
remanei <- "PX506"
briggsae <- colnames(ogs)[!(colnames(ogs) %in% nigoni) & colnames(ogs) != remanei & colnames(ogs)!="Orthogroup"]

og_qx_tran <- ogs %>% dplyr::select(Orthogroup,QX1410) %>%
  dplyr::filter(!grepl(",",QX1410)) %>%
  dplyr::filter(!is.na(QX1410)) 

briggsae_og_hdr <- lapply(briggsae[briggsae!="QX1410"], function(s) {
  
  hdr_dt <- hdrs %>%
    filter(STRAIN == s) %>%
    transmute(
      CHROM,
      start = minStart,
      end = maxEnd
    ) %>%
    data.table::as.data.table()
  
  tran_dt <- qxtran %>%
    filter(type == "mRNA") %>%
    mutate(
      transcript_id = stringr::str_match(attributes, "ID=transcript:([^;]+)")[, 2],
      og_id = paste0("transcript_", transcript_id),
      tran_len = end - start + 1
    ) %>%
    select(CHROM = seqid, start, end, tran_len, og_id) %>%
    data.table::as.data.table()
  
  data.table::setkey(hdr_dt, CHROM, start, end)
  data.table::setkey(tran_dt, CHROM, start, end)
  
  ov <- data.table::foverlaps(
    hdr_dt,
    tran_dt,
    nomatch =   0L
  )
  
  ov <- ov[
    ,
    overlap_bp := pmin(i.end, end) - pmax(i.start, start) + 1
  ][
    overlap_bp / tran_len > 0.50
  ]
  
  og_qx_tran %>%
    dplyr::mutate(
      og_id = stringr::str_remove(QX1410, "^QX1410\\|")
    ) %>%
    dplyr::semi_join(
      tibble::as_tibble(ov[, .(og_id)]),
      by = "og_id"
    ) %>%
    dplyr::transmute(
      OG = Orthogroup,
      STRAIN = s
    )
  
}) %>%
  bind_rows() %>%
  dplyr::distinct()

hdr_local <- hdrs %>%
  filter(STRAIN %in% briggsae) %>%
  dplyr::filter(STRAIN != "JU3200" & STRAIN != "QG4232") %>%
  transmute(
    CHROM,
    start = minStart,
    end = maxEnd
  ) %>%
  data.table::as.data.table()

get_strain <- function(x) sub("\\|.*$", "", x)

subset_aln <- function(aln, keep) {
  aln$seq <- aln$seq[keep]
  aln$nam <- aln$nam[keep]
  aln$nb  <- sum(keep)
  aln
}

read_and_filter_alignment <- function(aln_file) {
  aln <- seqinr::read.alignment(aln_file, format = "fasta")
  
  seqs <- toupper(aln$seq)
  
  occupancy <- vapply(seqs, function(x) {
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(chars %in% c("A", "C", "G", "T")) / length(chars)
  }, numeric(1))
  
  lenest <- vapply(seqs, function(x) {
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(chars %in% c("A", "C", "G", "T")) / 3
  }, numeric(1))
  
  keep <- occupancy >= 0.7 & lenest >= 100
  
  if (sum(keep) < 4) return(NULL)
  
  aln$nam <- aln$nam[keep]
  aln$seq <- aln$seq[keep]
  aln$nb <- length(aln$seq)
  
  aln
}


subset_alignment <- function(aln, keep) {
  aln$nam <- aln$nam[keep]
  aln$seq <- aln$seq[keep]
  aln$nb <- length(aln$seq)
  aln
}

extract_dist <- function(x, value_name) {
  m <- as.matrix(x)
  idx <- which(lower.tri(m), arr.ind = TRUE)
  
  data.frame(
    strain1 = rownames(m)[idx[, 1]],
    strain2 = colnames(m)[idx[, 2]],
    value = m[idx],
    stringsAsFactors = FALSE
  ) |>
    dplyr::rename(!!value_name := value)
}

sequinr_kaks_to_pair_df <- function(aln, og_id, comparison) {
  if (length(aln$seq) < 4) return(NULL)
  
  kk <- tryCatch(
    seqinr::kaks(aln),
    error = function(e) {
      warning("Skipping ", og_id, " ", comparison, ": kaks failed")
      return(NULL)
    }
  )
  
  if (is.null(kk)) return(NULL)
  
  ka_df  <- extract_dist(kk$ka,  "ka")
  ks_df  <- extract_dist(kk$ks,  "ks")
  vka_df <- extract_dist(kk$vka, "vka")
  vks_df <- extract_dist(kk$vks, "vks")
  
  out <- ka_df %>%
    dplyr::left_join(ks_df, by = c("strain1", "strain2")) %>%
    dplyr::left_join(vka_df, by = c("strain1", "strain2")) %>%
    dplyr::left_join(vks_df, by = c("strain1", "strain2")) %>%
    dplyr::mutate(
      kaks = dplyr::if_else(ks > 0, ka / ks, NA_real_),
      strain1 = get_strain(strain1),
      strain2 = get_strain(strain2),
      og_id = og_id,
      comparison = comparison
    ) 
  
  out
}


get_kaks_for_og_within_briggsae <- function(og_id) {
  
  suffixes <- c(
    ".no_cre.aligned.codons.gapped.fa",
    ".no_cn.aligned.codons.gapped.fa"
  )
  
  aln_file <- file.path(aln_dir, paste0(og_id, suffixes))
  aln_file <- aln_file[file.exists(aln_file)][1]
  
  if (is.na(aln_file)) {
    warning("Skipping ", og_id, ": no alignment found")
    return(NULL)
  }
  
  aln <- read_and_filter_alignment(aln_file)
  
  if (is.null(aln)) {
    warning("Skipping ", og_id, ": fewer than 4 sequences after filtering")
    return(NULL)
  }
  
  briggsae_keep <- get_strain(aln$nam) %in% briggsae
  
  if (sum(briggsae_keep) < 4) {
    warning("Skipping ", og_id, " within_briggsae: fewer than 4 briggsae sequences")
    return(NULL)
  }
  
  aln_briggsae <- subset_alignment(aln, briggsae_keep)
  
  sequinr_kaks_to_pair_df(
    aln_briggsae,
    og_id,
    comparison = "within_briggsae"
  ) %>%
    dplyr::mutate(pair_class = "within_briggsae")
}


get_kaks_for_og_nigoni_inter_species <- function(og_id) {
  
  aln_file <- file.path(
    aln_dir,
    paste0(og_id, ".no_cre.aligned.codons.gapped.fa")
  )
  
  if (!file.exists(aln_file)) {
    warning("Skipping ", og_id, ": no .no_cre alignment found")
    return(NULL)
  }
  
  aln <- read_and_filter_alignment(aln_file)
  
  if (is.null(aln)) {
    warning("Skipping ", og_id, ": fewer than 4 sequences after filtering")
    return(NULL)
  }
  
  full_df <- sequinr_kaks_to_pair_df(
    aln,
    og_id,
    comparison = "full_no_cre"
  )
  
  if (is.null(full_df)) return(NULL)
  
  full_df %>%
    dplyr::mutate(
      pair_class = dplyr::case_when(
        strain1 %in% briggsae & strain2 %in% briggsae ~ "within_briggsae_from_full",
        strain1 %in% briggsae | strain2 %in% briggsae ~ "inter_species",
        TRUE ~ "within_nigoni"
      )
    ) %>%
    dplyr::filter(pair_class %in% c("within_nigoni", "inter_species"))
}

hdr_within_briggsae <- bind_rows(
  lapply(hdr_groups, get_kaks_for_og_within_briggsae)
)

non_hdr_within_briggsae <- bind_rows(
  lapply(non_hdr_groups, get_kaks_for_og_within_briggsae)
)

hdr_all_spp <- bind_rows(
  lapply(hdr_groups, get_kaks_for_og_nigoni_inter_species)
)

non_hdr_all_spp <- bind_rows(
  lapply(non_hdr_groups, get_kaks_for_og_nigoni_inter_species)
)

hdr_pairwise_filtered <- hdr_within_briggsae%>%
  dplyr::mutate(class=ifelse(og_id %in% no_cre_files,"no_cre","no_cn")) %>%
  dplyr::mutate(
    sampleA = pmin(strain1, strain2),
    sampleB = pmax(strain1, strain2)
  ) %>%
  dplyr::filter(ks < 0.8) %>%
  dplyr::mutate(color=ifelse(sampleA == "QX1410" | sampleB == "QX1410","red","black"),
                comp=paste0(sampleA," vs. ",sampleB)) %>%
  dplyr::mutate(group="HDR")

non_hdr_pairwise_filtered <- non_hdr_within_briggsae%>%
  dplyr::mutate(class=ifelse(og_id %in% no_cre_files,"no_cre","no_cn")) %>%
  dplyr::mutate(
    sampleA = pmin(strain1, strain2),
    sampleB = pmax(strain1, strain2)
  ) %>%
  dplyr::filter(ks < 0.8) %>%
  dplyr::mutate(color=ifelse(sampleA == "QX1410" | sampleB == "QX1410","red","black"),
                comp=paste0(sampleA," vs. ",sampleB)) %>%
  dplyr::mutate(group="nHDR")

hdr_pairwise_filtered_all_spp <- hdr_all_spp%>%
  dplyr::mutate(class=ifelse(og_id %in% no_cre_files,"no_cre","no_cn")) %>%
  dplyr::mutate(
    sampleA = pmin(strain1, strain2),
    sampleB = pmax(strain1, strain2)
  ) %>%
  dplyr::filter(ks < 0.8) %>%
  dplyr::mutate(color=ifelse(sampleA == "QX1410" | sampleB == "QX1410","red","black"),
                comp=paste0(sampleA," vs. ",sampleB)) %>%
  dplyr::mutate(group="HDR")

non_hdr_pairwise_filtered_all_spp <- non_hdr_all_spp%>%
  dplyr::mutate(class=ifelse(og_id %in% no_cre_files,"no_cre","no_cn")) %>%
  dplyr::mutate(
    sampleA = pmin(strain1, strain2),
    sampleB = pmax(strain1, strain2)
  ) %>%
  dplyr::filter(ks < 0.8) %>%
  dplyr::mutate(color=ifelse(sampleA == "QX1410" | sampleB == "QX1410","red","black"),
                comp=paste0(sampleA," vs. ",sampleB)) %>%
  dplyr::mutate(group="nHDR")



hdr_pairwise_anno <- hdr_pairwise_filtered %>%
  dplyr::mutate(speciesA=ifelse(sampleA %in% nigoni,"nigoni","briggsae"),
                speciesB=ifelse(sampleB %in% nigoni,"nigoni","briggsae"))

non_hdr_pairwise_anno <- non_hdr_pairwise_filtered %>%
  dplyr::mutate(speciesA=ifelse(sampleA %in% nigoni,"nigoni","briggsae"),
                speciesB=ifelse(sampleB %in% nigoni,"nigoni","briggsae"))

hdr_pairwise_anno_all_spp <- hdr_pairwise_filtered_all_spp %>%
  dplyr::mutate(speciesA=ifelse(sampleA %in% nigoni,"nigoni","briggsae"),
                speciesB=ifelse(sampleB %in% nigoni,"nigoni","briggsae"))

non_hdr_pairwise_anno_all_spp <- non_hdr_pairwise_filtered_all_spp %>%
  dplyr::mutate(speciesA=ifelse(sampleA %in% nigoni,"nigoni","briggsae"),
                speciesB=ifelse(sampleB %in% nigoni,"nigoni","briggsae"))

all_pairwise_filtered <- rbind(hdr_pairwise_anno,non_hdr_pairwise_anno) %>%
  dplyr::mutate(species1=ifelse(sampleA %in% nigoni,"nigoni","briggsae"),
                species2=ifelse(sampleB %in% nigoni,"nigoni","briggsae")) %>%
  dplyr::mutate(comparison_check = case_when(
    species1 == species2 & species1 == "briggsae" ~ "within_briggsae",
    species1 == species2 & species1 == "nigoni"   ~ "within_nigoni",
    species1 != species2                          ~ "between_species",
    TRUE                                          ~ NA_character_
    )
  ) %>%
  dplyr::mutate(
    comparison_within = case_when(
      comparison_check == "within_briggsae" & ((sampleA=="JU3200" & sampleB!="QG4232") |(sampleB=="JU3200" & sampleA!="QG4232"))  ~ "Tropical_x_KD",
      comparison_check == "within_briggsae" & ((sampleA!="JU3200" & sampleB=="QG4232") |(sampleB!="JU3200" & sampleA=="QG4232"))  ~ "Tropical_x_AD",
      comparison_check == "within_briggsae" & sampleA!="JU3200" & sampleA!="QG4232" & sampleB!="QG4232" & sampleB!="JU3200"  ~ "Tropical_x_Tropical",
      comparison_check == "within_briggsae" & (sampleA=="JU3200" | sampleA=="QG4232") & (sampleB=="QG4232" | sampleB=="JU3200")  ~ "AD_x_KD",
      comparison_check == "within_nigoni" ~ "within_nigoni",
      comparison_check == "between_species" ~ "between_species",
      TRUE                                          ~ NA_character_
    )
  ) %>%
  dplyr::left_join(briggsae_og_hdr %>% dplyr::mutate(status1="HDR"), by=c("og_id"="OG","sampleA"="STRAIN")) %>%
  dplyr::left_join(briggsae_og_hdr %>% dplyr::mutate(status2="HDR"), by=c("og_id"="OG","sampleB"="STRAIN")) %>%
  dplyr::mutate(status1=ifelse(is.na(status1) & species1=="briggsae","non-HDR",status1)) %>%
  dplyr::mutate(status2=ifelse(is.na(status2) & species2=="briggsae","non-HDR",status2)) %>%
  dplyr::mutate(status1=ifelse(comparison_check=="between_species",NA,status1)) %>%
  dplyr::mutate(status2=ifelse(comparison_check=="between_species",NA,status2)) 


all_pairwise_filtered_all_spp <- rbind(hdr_pairwise_anno_all_spp,non_hdr_pairwise_anno_all_spp) %>%
  dplyr::mutate(species1=ifelse(sampleA %in% nigoni,"nigoni","briggsae"),
                species2=ifelse(sampleB %in% nigoni,"nigoni","briggsae")) %>%
  dplyr::mutate(comparison_check = case_when(
    species1 == species2 & species1 == "briggsae" ~ "within_briggsae",
    species1 == species2 & species1 == "nigoni"   ~ "within_nigoni",
    species1 != species2                          ~ "between_species",
    TRUE                                          ~ NA_character_
  )
  ) %>%
  dplyr::mutate(
    comparison_within = case_when(
      comparison_check == "within_briggsae" & ((sampleA=="JU3200" & sampleB!="QG4232") |(sampleB=="JU3200" & sampleA!="QG4232"))  ~ "Tropical_x_KD",
      comparison_check == "within_briggsae" & ((sampleA!="JU3200" & sampleB=="QG4232") |(sampleB!="JU3200" & sampleA=="QG4232"))  ~ "Tropical_x_AD",
      comparison_check == "within_briggsae" & sampleA!="JU3200" & sampleA!="QG4232" & sampleB!="QG4232" & sampleB!="JU3200"  ~ "Tropical_x_Tropical",
      comparison_check == "within_briggsae" & (sampleA=="JU3200" | sampleA=="QG4232") & (sampleB=="QG4232" | sampleB=="JU3200")  ~ "AD_x_KD",
      comparison_check == "within_nigoni" ~ "within_nigoni",
      comparison_check == "between_species" ~ "between_species",
      TRUE                                          ~ NA_character_
    )
  ) %>%
  dplyr::left_join(briggsae_og_hdr %>% dplyr::mutate(status1="HDR"), by=c("og_id"="OG","sampleA"="STRAIN")) %>%
  dplyr::left_join(briggsae_og_hdr %>% dplyr::mutate(status2="HDR"), by=c("og_id"="OG","sampleB"="STRAIN")) %>%
  dplyr::mutate(status1=ifelse(is.na(status1) & species1=="briggsae","non-HDR",status1)) %>%
  dplyr::mutate(status2=ifelse(is.na(status2) & species2=="briggsae","non-HDR",status2)) %>%
  dplyr::mutate(status1=ifelse(comparison_check=="between_species",NA,status1)) %>%
  dplyr::mutate(status2=ifelse(comparison_check=="between_species",NA,status2)) 

comp <- all_pairwise_filtered %>% 
  dplyr::filter(comparison_within=="Tropical_x_Tropical") %>%
  dplyr::mutate(comp_state=ifelse(status1==status2 & status1=="HDR","among_HDR",
                                  ifelse(status1==status2 & status1=="non-HDR","among_nHDR","across_hap"))) %>%
  dplyr::left_join(og_qx_tran %>% dplyr::mutate(tranname=gsub("QX1410\\|","",QX1410)) %>% dplyr::select(-QX1410),by=c("og_id"="Orthogroup")) %>%
  dplyr::left_join(qxtran %>%
                     tidyr::separate(attributes,into=c("ID","Parent","seqname","biotype","locus"),sep = ";") %>%
                     dplyr::mutate(ID=gsub(":","_",ID)) %>%
                     dplyr::mutate(ID=gsub("ID=","",ID)) %>%
                     dplyr::mutate(seqname=gsub("sequence_name=","",seqname)) %>%
                     dplyr::select(seqid,start,end,strand,ID,seqname), 
                   by=c("tranname"="ID")) %>%
  dplyr::group_by(og_id) %>%
  dplyr::mutate(
    group = if_else(
      any(status1 == "HDR" | status2 == "HDR"),
      "HDR",
      "non-HDR"
    )
  ) %>%
  dplyr::ungroup()


comp_plot <- comp %>%
  dplyr::group_by(sampleA, sampleB, comp_state) %>%
  dplyr::summarise(
    mean_ks = mean(ks, na.rm = TRUE),
    median_ks = median(ks, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(pairwise = paste(sampleA, sampleB, sep = " × ")) %>%
  tidyr::pivot_longer(
    cols = c(mean_ks, median_ks),
    names_to = "stat",
    values_to = "ks") %>%
  dplyr::mutate(stat = recode(stat,mean_ks = "Mean",median_ks="Median"))
# 

spcomp_plot <- all_pairwise_filtered_all_spp %>%
  dplyr::group_by(sampleA, sampleB, comparison_within) %>%
  dplyr::summarise(
    mean_ks = mean(ks, na.rm = TRUE),
    median_ks = median(ks, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(pairwise = paste(sampleA, sampleB, sep = " × ")) %>%
  tidyr::pivot_longer(
    cols = c(mean_ks, median_ks),
    names_to = "stat",
    values_to = "ks") %>%
  dplyr::mutate(stat = recode(stat,mean_ks = "Mean",median_ks="Median"))

spcomp_plot_by_group <- all_pairwise_filtered_all_spp %>%
  dplyr::group_by(sampleA, sampleB, comparison_within,group) %>%
  dplyr::summarise(
    mean_ks = mean(ks, na.rm = TRUE),
    median_ks = median(ks, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(pairwise = paste(sampleA, sampleB, sep = " × ")) %>%
  tidyr::pivot_longer(
    cols = c(mean_ks, median_ks),
    names_to = "stat",
    values_to = "ks") %>%
  dplyr::mutate(stat = recode(stat,mean_ks = "Mean",median_ks="Median"))

p2 <- spcomp_plot %>%
  dplyr::filter(stat == "Mean") %>%
  dplyr::mutate(
    comp_state_label = dplyr::recode(
      comparison_within,
      within_nigoni = "Among <i> C. nigoni</i>",
      between_species = paste(
        "Between <i>  C. briggsae</i><br>",
        "and <i>C. nigoni</i>",
        sep = " "
      )
    )
  ) %>%
  ggplot(aes(x = comp_state_label, y = ks)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) +
  ylab(expression(paste("Mean pairwise ", K[S]))) +
  ylim(0,0.22)


p3 <- spcomp_plot_by_group %>%
  dplyr::filter(stat == "Mean", comparison_within == "within_nigoni") %>%
  dplyr::mutate(
    comp_state_label = "Among <i>C. nigoni</i>"
  ) %>%
  ggplot(aes(x = group, y = ks)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~comp_state_label) +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(),  
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  scale_x_discrete(labels = c(HDR = "HDR",
                              nHDR = "non-HDR"  ))+
  ylab(expression(paste("Mean pairwise ", K[S]))) +
  ylim(0,NA)

ggsave(p3,filename = "../../figures/SF28_nigoni_hdr_divergence.png",device = "png",width = 7,height = 5,units = "in",bg = "white",dpi = 600)


comp2 <- all_pairwise_filtered %>% 
  dplyr::filter(comparison_within %in% c("Tropical_x_KD","Tropical_x_AD")) %>%
  dplyr::mutate(comp_state=ifelse(status1==status2 & status1=="HDR","among_HDR",
                                  ifelse(status1==status2 & status1=="non-HDR","among_nHDR","across_hap"))) %>%
  dplyr::left_join(og_qx_tran %>% dplyr::mutate(tranname=gsub("QX1410\\|","",QX1410)) %>% dplyr::select(-QX1410),by=c("og_id"="Orthogroup")) %>%
  dplyr::left_join(qxtran %>%
                     tidyr::separate(attributes,into=c("ID","Parent","seqname","biotype","locus"),sep = ";") %>%
                     dplyr::mutate(ID=gsub(":","_",ID)) %>%
                     dplyr::mutate(ID=gsub("ID=","",ID)) %>%
                     dplyr::mutate(seqname=gsub("sequence_name=","",seqname)) %>%
                     dplyr::select(seqid,start,end,strand,ID,seqname), 
                   by=c("tranname"="ID")) %>%
  dplyr::group_by(og_id) %>%
  #dplyr::mutate(mean_ks=mean(ks)) %>%
  dplyr::mutate(
    group = if_else(
      any(status1 == "HDR" | status2 == "HDR"),
      "HDR",
      "non-HDR"
    )
  ) %>%
  dplyr::ungroup()

comp2_plot <- comp2 %>%
  dplyr::group_by(sampleA, sampleB, comp_state,comparison_within) %>%
  dplyr::summarise(
    mean_ks = mean(ks, na.rm = TRUE),
    median_ks = median(ks, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(pairwise = paste(sampleA, sampleB, sep = " × ")) %>%
  tidyr::pivot_longer(
    cols = c(mean_ks, median_ks),
    names_to = "stat",
    values_to = "ks") %>%
  dplyr::mutate(stat = recode(stat,mean_ks = "Mean",median_ks="Median"))

p4 <- rbind(
  comp_plot %>%
    dplyr::filter(stat == "Mean") %>%
    dplyr::mutate(
      comparison_within = "Tropical_x_Tropical",
      meta = "Intra-population"
    ),
  comp2_plot %>%
    dplyr::filter(stat == "Mean") %>%
    dplyr::mutate(meta = "Inter-population")
) %>%
  dplyr::mutate(
    meta = factor(
      meta,
      levels = c("Intra-population", "Inter-population")
    ),
    comp_state_label = dplyr::recode(
      comp_state,
      "among_HDR" = "Among HDRs",
      "among_nHDR" = "Among non-HDRs",
      "across_hap" = "Between HDRs\nand non-HDRs"
    )
  ) %>%
  ggplot(aes(x = comp_state_label, y = ks)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(
    ~meta,
    nrow = 1,
    labeller = as_labeller(c(
      `Intra-population` = "Within Tropical RG",
      `Inter-population` = "Between RGs"
    ))
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_blank()
  ) +
  ylab(expression(paste("Mean pairwise ", K[S]))) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.005)
  ) +
  coord_cartesian(ylim = c(0, NA))


div_subplots <- cowplot::plot_grid(p4,p2,nrow=1,align = "h",axis = "tb",labels=c("b","c"),rel_widths = c(1,0.5))

saveRDS(div_subplots, file = "../../processed_data/tree_topology/aa_divergence/div_subplots.rds")


st14 <- rbind(comp_plot %>% 
                dplyr::mutate(comparison_within="Tropical_x_Tropical",
                              comparison="within_species") %>%
                dplyr::filter(stat=="Mean") %>%
                dplyr::rename(haplotypes=comp_state),
              comp2_plot %>%
                dplyr::filter(stat=="Mean") %>%
                dplyr::mutate(comparison="within_species") %>%
                dplyr::rename(haplotypes=comp_state),
              spcomp_plot_by_group %>%
                dplyr::filter(stat=="Mean") %>%
                dplyr::rename(comparison=comparison_within,
                              comparison_within=group) %>%
                dplyr::mutate(haplotypes=NA),
              spcomp_plot %>%
                dplyr::filter(stat=="Mean") %>%
                dplyr::rename(comparison=comparison_within) %>% 
                dplyr::mutate(comparison_within="All",
                              haplotypes=NA)) %>%
  dplyr::mutate(
    haplotypes_label = dplyr::case_when(
      haplotypes == "among_HDR" ~ "Among HDRs",
      haplotypes == "among_nHDR" ~ "Among non-HDRs",
      haplotypes == "across_hap" ~ "Between HDRs and non-HDRs",
      TRUE ~ haplotypes)) %>%
  dplyr::select(-stat,-pairwise,-haplotypes) %>%
  dplyr::select(sampleA,sampleB,comparison,comparison_within,haplotypes_label,ks) %>%
  dplyr::arrange(desc(comparison),desc(comparison_within),haplotypes_label,sampleA,sampleB) %>%
  dplyr::rename(`Mean Ks`=ks,Comparison=comparison,Comparison_within=comparison_within,`Sample A`=sampleA,`Sample B`=sampleB,Haplotypes=haplotypes_label) 

write.table(st14,"../../supplementary_data/ST14_mean_ks_est.tsv",sep = "\t",quote = F,row.names = F)
