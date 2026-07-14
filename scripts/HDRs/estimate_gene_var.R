library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)

hdrs <- readr::read_tsv("../../processed_data/HDRs/HDR_CB_allStrain_5kbclust_20260506.tsv") 
gff <- readr::read_tsv("../../processed_data/gene_diversity/c_briggsae.QX1410_20250929.csq.gff", col_names = c("seqid","source","type","start","end","score","strand","phase","attributes")) %>%
  dplyr::filter(seqid!="MtDNA")
ogs <- readr::read_tsv("../../processed_data/tree_topology/prot_TropicalCB_CN_CR_wDRG_39taxa/OrthoFinder/Results_CBCNCR_39_famsa/Orthogroups/Orthogroups.tsv")
ogs_colnames <- colnames(ogs)
colnames(ogs) <- sub("\\.longest\\.prot\\.otu$", "", ogs_colnames)
nigoni <- c("EG5268","JU2484","JU2617","JU4356","VSL2202","YR106","ZF1220","ECA2852","ECA2857","JU1418","JU1419","JU1420","JU1422","NIC2143","NIC2150","NIC2152","VX151","VX153")
remanei <- c("PX506")
adkd <- c("JU3200","QG4232")
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

tropical_ogs <- ogs %>%
  dplyr::select(-any_of(c(nigoni, remanei,adkd))) %>%
  dplyr::select(Orthogroup, QX1410, everything()) 
 

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

collapsed <- plyr::ldply(getRegFreq(hdrs %>% 
                                      dplyr::filter(source=="QX1410") %>%
                                      dplyr::group_split(CHROM)),data.frame) %>% 
  dplyr::mutate(divSize=maxEnd-minStart)

transcripts <- gff %>% dplyr::filter(type=="mRNA")

tx <- data.table::as.data.table(transcripts)
reg <- data.table::as.data.table(collapsed)

# foverlaps needs matching key column names
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

# set interval keys
data.table::setkey(reg, CHROM, reg_start, reg_end)

overlaps <- data.table::foverlaps(
  tx,
  reg,
  by.x = c("CHROM", "tx_start", "tx_end"),
  by.y = c("CHROM", "reg_start", "reg_end"),
  type = "any",
  nomatch = 0L
)

overlaps <- overlaps %>%
  dplyr::mutate(
    transcript_length = tx_end - tx_start + 1,
    overlap_start = pmax(tx_start, reg_start),
    overlap_end = pmin(tx_end, reg_end),
    overlap_bp = overlap_end - overlap_start + 1,
    pct_transcript_overlap = overlap_bp / transcript_length * 100
  )

overlap_tran <- overlaps %>% dplyr::filter(pct_transcript_overlap > 10) %>%
  dplyr::mutate(transcript=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(transcript=gsub("ID=","QX1410|",transcript)) %>%
  dplyr::mutate(transcript=gsub(":","_",transcript)) %>%
  dplyr::select(transcript,pct_transcript_overlap,CHROM,reg_start,reg_end,divSize)


tropical_ogs_overlap <- tropical_ogs %>% 
  tidyr::separate_rows(QX1410, sep = ",\\s*") %>%
  dplyr::left_join(overlap_tran,by=c("QX1410"="transcript"))

og_all_overlap <- tropical_ogs_overlap %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::filter(all(!is.na(pct_transcript_overlap))) %>%
  dplyr::ungroup()

og_any_overlap <- tropical_ogs_overlap %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::filter(any(!is.na(pct_transcript_overlap))) %>%
  dplyr::ungroup()

orthogroup_counts <- tropical_ogs %>%
  dplyr::mutate(dplyr::across(-Orthogroup,~ dplyr::case_when(is.na(.) ~ 0,TRUE ~ stringr::str_count(., ",") + 1))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    classification = dplyr::case_when(
      all(c_across(-c(Orthogroup)) == 0) ~ "NCB",
      all(c_across(-c(Orthogroup)) == 1) ~ "SIC",
      length(unique(c_across(-c(Orthogroup)))) == 1 ~ "STC",
      TRUE ~ "CNV"
    )
  ) %>%
  dplyr::ungroup()
                              
domains <- data.table::as.data.table(domains_wide)

# Build ARM intervals: left arm and right arm
arms <- data.table::rbindlist(
  list(
    domains[, .(
      CHROM,
      domain = "LEFT_ARM",
      arm_start = left_start,
      arm_end = left_end
    )],
    domains[, .(
      CHROM,
      domain = "RIGHT_ARM",
      arm_start = right_start,
      arm_end = right_end
    )]
  )
)

data.table::setkey(arms, CHROM, arm_start, arm_end)

arm_overlaps <- data.table::foverlaps(
  tx,
  arms,
  by.x = c("CHROM", "tx_start", "tx_end"),
  by.y = c("CHROM", "arm_start", "arm_end"),
  type = "any",
  nomatch = 0L
)

arm_overlaps <- arm_overlaps %>%
  dplyr::mutate(
    transcript_length = tx_end - tx_start + 1,
    overlap_start = pmax(tx_start, arm_start),
    overlap_end = pmin(tx_end, arm_end),
    overlap_bp = overlap_end - overlap_start + 1,
    pct_arm_overlap = overlap_bp / transcript_length * 100
  ) %>%
  dplyr::filter(pct_arm_overlap > 10) %>%
  dplyr::distinct(CHROM, tx_start, tx_end, .keep_all = TRUE) %>%
  dplyr::select(CHROM, tx_start, tx_end, domain, pct_arm_overlap)

tx_classified <- tx %>%
  dplyr::left_join(
    arm_overlaps,
    by = c("CHROM", "tx_start", "tx_end")
  ) %>%
  dplyr::left_join(
    domains,
    by = "CHROM"
  ) %>%
  dplyr::mutate(
    region = dplyr::case_when(
      !is.na(pct_arm_overlap) ~ "ARM",
      tx_end < left_start | tx_start > right_end ~ "TIP",
      tx_start > left_end & tx_end < right_start ~ "CENTER",
      TRUE ~ "BOUNDARY"
    )
  ) %>%
  dplyr::mutate(transcript=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(transcript=gsub("ID=","QX1410|",transcript)) %>%
  dplyr::mutate(transcript=gsub(":","_",transcript)) %>%
  dplyr::select(transcript,pct_arm_overlap,CHROM,region)

#confirm classification visually
#ggplot(tx_classified) +geom_point(aes(x=tx_start, y=1, color=region)) + facet_wrap(~CHROM,ncol=1)
tropical_ogs_overlap_dom <- tropical_ogs %>% 
  tidyr::separate_rows(QX1410, sep = ",\\s*") %>%
  dplyr::left_join(tx_classified,by=c("QX1410"="transcript"))

tropical_og_overlap_dom_class <- tropical_ogs_overlap_dom %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::mutate(
    domain_class = dplyr::case_when(
      all(region == "ARM", na.rm = FALSE) ~ "FULL_ARM",
      all(region == "TIP", na.rm = FALSE) ~ "FULL_TIP",
      all(region == "CENTER", na.rm = FALSE) ~ "FULL_CENTER",
      TRUE ~ "MIXED"
    )
  ) %>%
  dplyr::ungroup()

tropical_og_overlap_dom_summary <- tropical_og_overlap_dom_class %>%
  dplyr::distinct(Orthogroup, domain_class)

all_present <- orthogroup_counts %>% 
  dplyr::filter(classification != "NCB") %>%
  dplyr::mutate(overlap=ifelse(Orthogroup %in% og_all_overlap$Orthogroup,"FULL",
                               ifelse(Orthogroup %in% og_any_overlap$Orthogroup,"PARTIAL","NONE"))) %>%
  dplyr::left_join(tropical_og_overlap_dom_summary,by="Orthogroup") %>%
  dplyr::rename(hdr_overlap = overlap, domain_overlap=domain_class)


hdr_og_list <- unique((all_present %>% dplyr::filter(hdr_overlap =="FULL"))$Orthogroup)
n_hdr_og_list <- unique((all_present %>% dplyr::filter(hdr_overlap =="NONE"))$Orthogroup)

writeLines(hdr_og_list, "../../processed_data/gene_diversity/hdr_og_list.txt")
writeLines(n_hdr_og_list, "../../processed_data/gene_diversity/non_hdr_og_list.txt")

enrichment_test <- function(
    df,
    class_value,
    overlap_values,
    overlap_col = "hdr_overlap"
) {
  overlap_vec <- df[[overlap_col]]
  
  N <- nrow(df)
  O <- sum(overlap_vec %in% overlap_values)
  C <- sum(df$classification == class_value)
  K <- sum(
    df$classification == class_value &
      overlap_vec %in% overlap_values
  )
  
  mat <- matrix(
    c(
      K,
      C - K,
      O - K,
      N - O - C + K
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Class = c(class_value, paste0("Not_", class_value)),
      Overlap = c(
        paste(overlap_values, collapse = "_or_"),
        paste0("Not_", paste(overlap_values, collapse = "_or_"))
      )
    )
  )
  
  fisher_out <- fisher.test(mat, alternative = "two.sided")
  
  tibble::tibble(
    class = class_value,
    overlap_set = paste(overlap_values, collapse = "+"),
    N = N,
    O = O,
    C = C,
    K = K,
    fold_enrichment = (K / C) / (O / N),
    odds_ratio = unname(fisher_out$estimate),
    p_value = fisher_out$p.value,
    conf.low = fisher_out$conf.int[1],
    conf.high = fisher_out$conf.int[2],
    direction = ifelse(fold_enrichment > 1, "enriched", "depleted")
  )
}

classes <- c("CNV", "STC", "SIC")

overlap_sets <- list(
  FULL = "FULL",
  FULL_PARTIAL = c("FULL", "PARTIAL")
)

background_sets <- list(
  ALL_ORTHOGROUPS = all_present,
  FULL_ARM_ONLY = all_present %>%
    dplyr::filter(domain_overlap == "FULL_ARM")
)

enrichment_results <- purrr::map_dfr(
  names(background_sets),
  function(background_name) {
    purrr::map_dfr(
      classes,
      function(class_name) {
        purrr::map_dfr(
          names(overlap_sets),
          function(overlap_name) {
            enrichment_test(
              df = background_sets[[background_name]],
              class_value = class_name,
              overlap_values = overlap_sets[[overlap_name]],
              overlap_col = "hdr_overlap"
            ) %>%
              dplyr::mutate(
                background = background_name,
                overlap_group = overlap_name,
                .before = overlap_set
              )
          }
        )
      }
    )
  }
)

enrichment_results_adj <- enrichment_results %>%
  dplyr::mutate(
    p_adj_bh = p.adjust(p_value, method = "BH"),
    p_adj_bonferroni = p.adjust(p_value, method = "bonferroni")
  ) %>%
  dplyr::mutate(gene_class=ifelse(class=="STC","stable_copy_ortholog",ifelse(class=="CNV","variable_copy_ortholog","single_copy_ortholog"))) %>%
  dplyr::select(-class, -overlap_group) %>%
  dplyr::rename(total_orthogroups = N,
                hdr_overlapping_orthogroups = O,
                gene_class_orthogroups = C,
                overlapping_gene_class_orthogroups = K) %>%
  dplyr::select(ortholog_class=gene_class, everything()) %>%
  dplyr::arrange(ortholog_class)

write.table(enrichment_results_adj, file = "../../supplementary_data/ST9_ortholog_copy_enrichment_tests.tsv",row.names = F,quote = F,sep = "\t")
