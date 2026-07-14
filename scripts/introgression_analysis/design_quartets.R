library(ape)
library(ggplot2)
library(ggtree)
library(cowplot)
library(dplyr)
library(tidyr)

hdrs <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_allStrain_5kbclust_20260506.tsv") 
gff <- readr::read_tsv("../../processed_data/gene_diversity/c_briggsae.QX1410_20250929.csq.gff", col_names = c("seqid","source","type","start","end","score","strand","phase","attributes")) %>%
  dplyr::filter(seqid!="MtDNA")
transcripts <- gff %>% dplyr::filter(type=="mRNA")
ogs <- readr::read_tsv("../../processed_data/tree_topology/prot_TropicalCB_CN_CR_wDRG_39taxa/OrthoFinder/Results_CBCNCR_39_famsa/Orthogroups/Orthogroups.tsv")
ogs_colnames <- colnames(ogs)
colnames(ogs) <- sub("\\.longest\\.prot\\.otu$", "", ogs_colnames)
nigoni <- c("EG5268","JU2484","JU2617","YR106","ZF1220","ECA2852","ECA2857","JU1418","JU1419","JU1420","JU1422","NIC2143","NIC2150","NIC2152","VX151","VX153") 
remanei <- "PX506"
briggsae <- colnames(ogs)[!(colnames(ogs) %in% nigoni) & colnames(ogs) != remanei & colnames(ogs)!="Orthogroup"]

sample_cols <- c(nigoni, briggsae)
all_samp <- c(nigoni,remanei,briggsae)

og_counts_all <- ogs %>%
  dplyr::mutate(QX1410_original = QX1410) %>%
  dplyr::select(Orthogroup, QX1410_original, all_of(all_samp)) %>%
  dplyr::mutate(
    dplyr::across(
      all_of(all_samp),
      ~ dplyr::if_else(
        is.na(.),
        NA_integer_,
        stringr::str_count(., ",") + 1L
      )
    )
  )

#####EXTRACT COUNTS########
###########################
# all_na <- function(x) all(is.na(x))
# any_present <- function(x) any(!is.na(x))
# all_one <- function(x) all(x[!is.na(x)] == 1)
# full_one <- function(x) all(!is.na(x) & x == 1)
# 
# summary_df <- og_counts_all %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(
#     ## species-specific presence
#     nig_present = any_present(c_across(all_of(nigoni))),
#     brig_present = any_present(c_across(all_of(briggsae))),
#     rem_present = !is.na(c_across(all_of(remanei))),
#     
#     ## species-specific absence
#     nig_absent = all_na(c_across(all_of(nigoni))),
#     brig_absent = all_na(c_across(all_of(briggsae))),
#     rem_absent = is.na(c_across(all_of(remanei))),
#     
#     ## full occupancy within species
#     nig_full_one = full_one(c_across(all_of(nigoni))),
#     brig_full_one = full_one(c_across(all_of(briggsae))),
#     rem_one = c_across(all_of(remanei)) == 1
#   ) %>%
#   dplyr::ungroup()
# 
# unique_nigoni <- summary_df %>%
#   filter(nig_present, brig_absent, rem_absent)
# 
# unique_briggsae <- summary_df %>%
#   filter(brig_present, nig_absent, rem_absent)
# 
# unique_remanei <- summary_df %>%
#   filter(rem_present, nig_absent, brig_absent)
# 
# singlecopy_nigoni <- og_counts_all %>%
#   rowwise() %>%
#   filter(all_one(c_across(all_of(nigoni))))
# 
# singlecopy_briggsae <- og_counts_all %>%
#   rowwise() %>%
#   filter(all_one(c_across(all_of(briggsae))))
# 
# singlecopy_remanei <- og_counts_all %>%
#   filter(.data[[remanei]] == 1)
# 
# shared_three <- summary_df %>%
#   filter(nig_present, brig_present, rem_present)
# 
# shared_single <- og_counts_all %>%
#   rowwise() %>%
#   filter(
#     any_present(c_across(all_of(nigoni))),
#     any_present(c_across(all_of(briggsae))),
#     !is.na(c_across(all_of(remanei))),
#     all_one(c_across(all_of(nigoni))),
#     all_one(c_across(all_of(briggsae))),
#     c_across(all_of(remanei)) == 1
#   )
# 
# shared_single_full <- og_counts_all %>%
#   rowwise() %>%
#   filter(
#     full_one(c_across(all_of(nigoni))),
#     full_one(c_across(all_of(briggsae))),
#     c_across(all_of(remanei)) == 1
#   )
# 
# summary_tbl <- tibble(
#   category = c(
#     "Unique to C. nigoni",
#     "Unique to C. briggsae",
#     "Unique to C. remanei",
#     "Single-copy C. nigoni",
#     "Single-copy C. briggsae",
#     "Single-copy C. remanei",
#     "Shared across all species",
#     "Shared single-copy",
#     "Shared single-copy (full occupancy)"
#   ),
#   n_orthogroups = c(
#     nrow(unique_nigoni),
#     nrow(unique_briggsae),
#     nrow(unique_remanei),
#     nrow(singlecopy_nigoni),
#     nrow(singlecopy_briggsae),
#     nrow(singlecopy_remanei),
#     nrow(shared_three),
#     nrow(shared_single),
#     nrow(shared_single_full)
#   )
# )
#############################

all_sc_orthos <- ogs %>%
  dplyr::mutate(QX1410_original = QX1410) %>%
  dplyr::select(Orthogroup, QX1410_original, all_of(sample_cols)) %>%
  dplyr::mutate(
    dplyr::across(
      all_of(sample_cols),
      ~ dplyr::if_else(
        is.na(.),
        NA_integer_,
        stringr::str_count(., ",") + 1L
      )
    )
  ) %>%
  dplyr::filter(!dplyr::if_any(all_of(briggsae), ~ !is.na(.) & . > 1)) %>%
  dplyr::filter(rowSums(!is.na(dplyr::pick(all_of(briggsae)))) >= 3) %>%
  dplyr::filter(!is.na(QX1410_original))


all_sc_orthos_wcn <- all_sc_orthos %>%
  dplyr::filter(!dplyr::if_any(all_of(sample_cols), ~ !is.na(.) & . > 1)) %>%
  dplyr::pull(Orthogroup)


focal_cb_strains <- hdrs %>% 
  dplyr::filter(STRAIN %in% briggsae) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(span=sum(divSize)) %>%
  dplyr::filter(span==max(span))

all_focal_cb <- c("JU3200","JU3237","QG4232")
all_potential_scorthos_wcncr <- og_counts_all %>%
  dplyr::filter(
    rowSums(!is.na(pick(all_of(briggsae)))) >= 3,
    if_any(all_of(all_focal_cb), ~ . == 1),
    if_any(all_of(nigoni), ~ . == 1),
    if_any(all_of(remanei), ~ . == 1),
    QX1410==1
  ) #%>%
  dplyr::select(Orthogroup,QX1410_original,QX1410,PX506,all_of(c(all_focal_cb,nigoni))) %>%
  dplyr::pull(Orthogroup)

trim_aln <- all_samp[!(all_samp %in% c(nigoni,remanei,all_focal_cb,"QX1410"))]
writeLines(all_potential_scorthos_wcncr, "../../processed_data/tree_topology/phylo_py_pruner/all_pot_paired_file_groups.tsv")
writeLines(trim_aln, "../../processed_data/tree_topology/phylo_py_pruner/samples_to_trim_quartet.txt")

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
        
        retain <- temp %>% dplyr::filter(check == F & dplyr::coalesce(dplyr::lag(check), F) == F)
        
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

focal_cb_strain <- "QG4232"
focal_hdrs <- hdrs %>% dplyr::filter(STRAIN==focal_cb_strain)
tropical_hdrs<- plyr::ldply(getRegFreq(hdrs %>% 
                                         dplyr::filter(STRAIN %in% briggsae) %>%
                                         dplyr::group_split(CHROM)),data.frame) %>% 
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::select(-STRAIN)
                                                    

tx <- data.table::as.data.table(transcripts)
reg <- data.table::as.data.table(focal_hdrs)
reg_trop <- data.table::as.data.table(tropical_hdrs)

data.table::setnames(
  tx,
  c("seqid", "start", "end"),
  c("CHROM", "tx_start", "tx_end")
)

data.table::setnames(
  reg,
  c("minStart", "maxEnd"),
  c("reg_start", "reg_end")
)

data.table::setnames(
  reg_trop,
  c("minStart", "maxEnd"),
  c("reg_start", "reg_end")
)

data.table::setkey(reg, CHROM, reg_start, reg_end)
data.table::setkey(reg_trop, CHROM, reg_start, reg_end)

overlaps <- data.table::foverlaps(
  tx,
  reg,
  by.x = c("CHROM", "tx_start", "tx_end"),
  by.y = c("CHROM", "reg_start", "reg_end"),
  type = "any",
  nomatch = 0L
)

trop_overlaps <- data.table::foverlaps(
  tx,
  reg_trop,
  by.x = c("CHROM", "tx_start", "tx_end"),
  by.y = c("CHROM", "reg_start", "reg_end"),
  type = "any",
  nomatch = 0L
)

overlaps_wstats <- overlaps %>%
  dplyr::mutate(
    transcript_length = tx_end - tx_start + 1,
    overlap_start = pmax(tx_start, reg_start),
    overlap_end = pmin(tx_end, reg_end),
    overlap_bp = overlap_end - overlap_start + 1,
    pct_transcript_overlap = overlap_bp / transcript_length * 100
  )

trop_overlaps_wstats <- trop_overlaps %>%
  dplyr::mutate(
    transcript_length = tx_end - tx_start + 1,
    overlap_start = pmax(tx_start, reg_start),
    overlap_end = pmin(tx_end, reg_end),
    overlap_bp = overlap_end - overlap_start + 1,
    pct_transcript_overlap = overlap_bp / transcript_length * 100
  )

overlap_tran <- overlaps_wstats %>% dplyr::filter(pct_transcript_overlap > 10) %>%
  dplyr::mutate(transcript=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(transcript=gsub("ID=","QX1410|",transcript)) %>%
  dplyr::mutate(transcript=gsub(":","_",transcript)) %>%
  dplyr::select(transcript,pct_transcript_overlap,CHROM,reg_start,reg_end,divSize)

trop_overlap_tran <- trop_overlaps_wstats %>% dplyr::filter(pct_transcript_overlap > 10) %>%
  dplyr::mutate(transcript=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(transcript=gsub("ID=","QX1410|",transcript)) %>%
  dplyr::mutate(transcript=gsub(":","_",transcript)) %>%
  dplyr::select(transcript,pct_transcript_overlap,CHROM,reg_start,reg_end,divSize)

ogs_tran_merge <- ogs%>%
  dplyr::select(Orthogroup, QX1410, everything()) %>%
  tidyr::separate_rows(QX1410, sep = ",\\s*") %>%
  dplyr::left_join(overlap_tran,by=c("QX1410"="transcript"))

trop_ogs_tran_merge <- all_sc_orthos %>%
  dplyr::select(Orthogroup, QX1410_original) %>%
  dplyr::rename(QX1410=QX1410_original) %>%
  tidyr::separate_rows(QX1410, sep = ",\\s*") %>%
  dplyr::left_join(trop_overlap_tran,by=c("QX1410"="transcript"))

og_all_overlap <- ogs_tran_merge %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::filter(all(!is.na(pct_transcript_overlap))) %>%
  dplyr::ungroup()

trop_og_all_overlap <- trop_ogs_tran_merge %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::filter(all(!is.na(pct_transcript_overlap))) %>%
  dplyr::ungroup()

og_any_overlap <- ogs_tran_merge %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::filter(any(!is.na(pct_transcript_overlap))) %>%
  dplyr::ungroup()

trop_og_any_overlap <- trop_ogs_tran_merge %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::filter(any(!is.na(pct_transcript_overlap))) %>%
  dplyr::ungroup()

ogs_classified <- ogs %>% 
  dplyr::mutate(overlap=ifelse(Orthogroup %in% og_all_overlap$Orthogroup,"FULL",
                               ifelse(Orthogroup %in% og_any_overlap$Orthogroup,"PARTIAL","NONE"))) %>%
  dplyr::rename(hdr_overlap = overlap) 

trop_ogs_classified <- ogs %>% 
  dplyr::filter(Orthogroup %in% all_sc_orthos$Orthogroup) %>%
  dplyr::mutate(overlap=ifelse(Orthogroup %in% trop_og_all_overlap$Orthogroup,"FULL",
                               ifelse(Orthogroup %in% trop_og_any_overlap$Orthogroup,"PARTIAL","NONE"))) %>%
  dplyr::rename(hdr_overlap = overlap) 

hdr_ogs <- ogs_classified %>% 
  dplyr::filter(hdr_overlap=="FULL") %>%
  dplyr::filter(Orthogroup %in% all_potential_scorthos_wcncr) %>%
  dplyr::pull(Orthogroup)

non_hdr_ogs <- ogs_classified %>% 
  dplyr::filter(hdr_overlap=="NONE") %>%
  dplyr::filter(Orthogroup %in% all_potential_scorthos_wcncr) %>%
  dplyr::pull(Orthogroup)

writeLines(hdr_ogs, paste0("../../processed_data/tree_topology/phylo_py_pruner/hdr_ogs_",focal_cb_strain,".txt"))
writeLines(non_hdr_ogs, paste0("../../processed_data/tree_topology/phylo_py_pruner/non_hdr_ogs_",focal_cb_strain,".txt"))

trop_hdr_ogs <- trop_ogs_classified %>% 
  dplyr::filter(hdr_overlap=="FULL") %>%
  dplyr::pull(Orthogroup)

trop_non_hdr_ogs <- trop_ogs_classified %>% 
  dplyr::filter(hdr_overlap=="NONE") %>%
  dplyr::pull(Orthogroup)



all_trop_sc_cbonly <- c(trop_hdr_ogs,trop_non_hdr_ogs)
all_trop_sc_cbcn <- all_trop_sc_cbonly[all_trop_sc_cbonly %in% all_sc_orthos_wcn]

writeLines(trop_hdr_ogs, paste0("../../processed_data/tree_topology/aa_divergence/hdr_ogs_all_tropical.txt"))
writeLines(trop_non_hdr_ogs, paste0("../../processed_data/tree_topology/aa_divergence/non_hdr_ogs_all_tropical.txt"))
writeLines(all_trop_sc_cbonly, paste0("../../processed_data/tree_topology/aa_divergence/all_tropical_classified_ogs.txt"))
writeLines(all_trop_sc_cbcn, paste0("../../processed_data/tree_topology/aa_divergence/all_tropical_classified_scall.txt"))


ref <- "QX1410" 
focal_cb_strain <- "JU3237"
quartet_df <- do.call(rbind, lapply(nigoni, function(cnx) {
  data.frame(
    CNX = cnx,
    topology = c("O_CB1CB2_CNX", "O_CB1CNX_CB2", "O_CB2CNX_CB1"),
    quartet = c(
      sprintf("(%s,(%s,%s),%s)", remanei, ref, focal_cb_strain, cnx),
      sprintf("(%s,(%s,%s),%s)", remanei, ref, cnx, focal_cb_strain),
      sprintf("(%s,(%s,%s),%s)", remanei, focal_cb_strain, cnx, ref)
    ),
    samples = rep(
      paste(c(remanei, ref, focal_cb_strain, cnx), collapse = ","),
      3
    ),
    stringsAsFactors = FALSE
  )
}))


quartet_df_enum <- quartet_df %>%
  dplyr::group_by(samples) %>%
  dplyr::mutate(
    quartet_id = sprintf("quartet_%02d", dplyr::cur_group_id())
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(quartet_id) %>%
  dplyr::mutate(
    exclude = purrr::map_chr(
      samples,
      ~ paste(
        setdiff(all_samp, strsplit(.x, ",", fixed = TRUE)[[1]]),
        collapse = ","
      )
    )
  )

outdir <- "../../processed_data/tree_topology/mast/quartet_topologies"

quartet_df_enum %>%
  dplyr::group_by(quartet_id) %>%
  dplyr::group_walk(~ {
    outfile <- file.path(
      outdir,
      paste0(
        .y$quartet_id,
        "_",
        dplyr::first(.x$CNX),
        "_",
        focal_cb_strain,
        "_topologies.tre"
      )
    )
    
    writeLines(
      paste0(.x$quartet, ";"),
      outfile
    )
  })

quartet_exclude <- quartet_df_enum %>%
  dplyr::distinct(quartet_id, exclude) %>%
  dplyr::mutate(
    exclude_otus = purrr::map_chr(
      strsplit(exclude, ",", fixed = TRUE),
      ~ paste0(.x, "_longest_prot_otu_", .x, collapse = " ")
    )
  )

write.table( quartet_exclude,paste0("../../processed_data/tree_topology/phylo_py_pruner/",focal_cb_strain,"_quartet_exclude.tsv"),quote = F,sep = "\t",row.names = F)
write.table( quartet_df_enum %>% dplyr::select(-exclude),paste0("../../processed_data/tree_topology/mast/quartet_topologies/summ_quartet_topos/",focal_cb_strain,"_quartets_summary.tsv"),quote = F,sep = "\t",row.names = F)