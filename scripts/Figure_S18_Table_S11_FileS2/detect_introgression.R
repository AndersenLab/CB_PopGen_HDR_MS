#library('TreeDist')
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridGraphics)
library(tibble)
library(ggrepel)
library(ggtree)
library(phangorn)
library(data.table)
library(purrr)
library(circlize)
library(ComplexHeatmap)
library(stringr)

process_sites <- function(path, label) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    dplyr::filter(!grepl("MtDNA", POS)) %>%
    # compute missingness across all genotype columns (everything except the POS string col)
    dplyr::mutate(
      n_missing = rowSums(as.data.frame(dplyr::across(-POS, ~ is.na(.x)))),
      n_samples = ncol(.) - 1L,
      pct_missing = n_missing / n_samples
    ) %>%
    tidyr::separate(POS, into = c("CHROM", "POS"), sep = ":", remove = FALSE) %>%
    dplyr::mutate(POS = readr::parse_number(POS)) %>%
    dplyr::select(CHROM, POS, n_samples, n_missing, pct_missing) %>%
    dplyr::mutate(PL = label)
}

orthogroups <- readr::read_tsv("../processed_data/Orthogroups/Orthogroups.tsv")

#file_list <- list.files(path = "/vast/eande106/projects/Nicolas/CBCN_comp/tree/gene_level_20251006/core_prot_CBCN_fasttree/OrthoFinder/Results_CBCN_famsa_asn_ft/", pattern = "\\.txt$", full.names = TRUE)

tree_lines <- readLines("../processed_data/Resolved_Gene_Trees.txt")
trees_raw <- grep("^OG[0-9]+:", tree_lines, value = TRUE)
tree_texts <- sub("^OG[0-9]+:\\s*", "", trees_raw)
trees <- lapply(tree_texts, function(x) read.tree(text = x))
og_ids <- sub(":.*", "", grep("^OG[0-9]+:", tree_lines, value = TRUE))
names(trees) <- og_ids
trees <- lapply(trees, function(tr) {
  tr$tip.label <- gsub("_longest_prot", "", tr$tip.label)
  tr$tip.label <- gsub("_Transcript", "", tr$tip.label)
  tr$tip.label <- gsub("_transcript", "", tr$tip.label)
  tr
})
#glimpse(trees)

consensus <- ape::read.tree(file="../processed_data/SpeciesTree_rooted_node_labels.txt")
consensus$tip.label <- gsub(".longest.prot","",consensus$tip.label)

all_str <- consensus$tip.label
nigoni <- c("EG5268","JU2484","JU2617","JU4356","VSL2202","YR106","ZF1220","ECA2852","ECA2857","JU1418","JU1419","JU1420","JU1422","NIC2143","NIC2150","NIC2152","VX151","VX153")
briggsae <- all_str[!(all_str %in% nigoni)]

c1 <- colnames(orthogroups) 
c2 <- gsub("\\.longest\\.prot","",c1)
colnames(orthogroups) <- c2


og <- orthogroups

# Keep only strains that actually exist as columns
present_nigoni   <- intersect(nigoni, names(og))
present_briggsae <- intersect(briggsae, names(og))

og_classified <- og %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    # store per-row vectors as list-cols
    vals_nig = list(c_across(any_of(present_nigoni))),
    vals_bri = list(c_across(any_of(present_briggsae))),
    
    # unlist versions
    u_nig = list(unlist(vals_nig)),
    u_bri = list(unlist(vals_bri)),
    
    nigoni_has   = any(!is.na(u_nig)),
    briggsae_has = any(!is.na(u_bri)),
    
    nigoni_multi   = any(str_detect(u_nig, ","), na.rm = TRUE),
    briggsae_multi = any(str_detect(u_bri, ","), na.rm = TRUE),
    
    both_vals   = list(c(u_nig, u_bri)),
    
    # all cells (that are non-NA) have no commas
    all_single  = {
      x <- unlist(both_vals)
      all(is.na(x) | !str_detect(x, ","))
    },
    # every cell present (no NAs) AND we actually had columns
    all_present = {
      x <- unlist(both_vals)
      length(x) > 0 && all(!is.na(x))
    }
  ) %>%
  dplyr::mutate(
    group_type = case_when(
      nigoni_has   & !briggsae_has ~ "nigoni_unique",
      briggsae_has & !nigoni_has   ~ "briggsae_unique",
      nigoni_has   &  briggsae_has ~ "shared",
      TRUE ~ "none"
    ),
    # shared + no commas in any non-NA cell for both species
    shared_single_any = group_type == "shared" & !nigoni_multi & !briggsae_multi,
    # shared + all present (no NAs) + all singletons (no commas)
    shared_single_all = group_type == "shared" & all_present & all_single
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-u_nig, -u_bri)

# Subsets
unique_nigoni     <- filter(og_classified, group_type == "nigoni_unique")
unique_briggsae   <- filter(og_classified, group_type == "briggsae_unique")
shared_groups     <- filter(og_classified, group_type == "shared")
shared_single_any <- filter(og_classified, shared_single_any)
shared_single_all <- filter(og_classified, shared_single_all)

# Summaries
table(og_classified$group_type)
sum(og_classified$shared_single_any)
sum(og_classified$shared_single_all)
############

sample_from_label <- function(lbl) sub("_.*$", "", lbl)

all_samples <- c(nigoni, briggsae)

label_to_sample <- function(lbl) {
  smp <- sample_from_label(lbl)
  if (smp %in% all_samples) smp else NA_character_
}

label_to_species <- function(lbl) {
  smp <- label_to_sample(lbl)
  if (is.na(smp)) return(NA_character_)
  if (smp %in% nigoni)   return("nigoni")
  if (smp %in% briggsae) return("briggsae")
  NA_character_
}

# data.frame(
#   label   = example_labels,
#   sample  = vapply(example_labels, label_to_sample, character(1)),
#   species = vapply(example_labels, label_to_species, character(1)),
#   stringsAsFactors = FALSE
# )

species_of <- function(tips) {
  out <- vapply(tips, label_to_species, character(1))
  names(out) <- tips
  out
}

is_mono <- function(tr, which_species) {
  sp_map  <- species_of(tr$tip.label)
  tips_in <- names(sp_map)[sp_map == which_species]
  if (length(tips_in) < 2) return(NA)   # not enough to test
  is.monophyletic(tr, tips_in)
}

# Find *minimal* mixed clades (smallest subtrees containing both species)
mixed_clades <- function(tr) {
  sp_map <- species_of(tr$tip.label)
  if (sum(!is.na(sp_map)) < 4) return(list())
  
  mixes <- list()
  nodes <- (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode)
  k <- 0
  for (nd in nodes) {
    desc <- Descendants(tr, nd, type = "tips")[[1]]
    spd  <- na.omit(sp_map[tr$tip.label[desc]])
    if (length(unique(spd)) > 1) {
      # minimality: no child also mixed
      mixed_child <- FALSE
      for (ch in Children(tr, nd)) {
        dch  <- Descendants(tr, ch, type = "tips")[[1]]
        spch <- na.omit(sp_map[tr$tip.label[dch]])
        if (length(unique(spch)) > 1) { mixed_child <- TRUE; break }
      }
      if (!mixed_child) {
        k <- k + 1
        mixes[[k]] <- list(
          node = nd,
          tips = tr$tip.label[desc],
          species_counts = table(spd),
          samples = unique(sample_from_label(tr$tip.label[desc])) # which samples involved
        )
      }
    }
  }
  mixes
}

is_single_copy_tree <- function(tr) {
  if (!inherits(tr, "phylo")) return(NA)
  ids <- vapply(tr$tip.label, sample_from_label, character(1))
  ids <- ids[ids %in% all_samples]              # only consider your samples
  if (length(ids) == 0) return(NA)              # nothing to evaluate
  !any(duplicated(ids))                         # TRUE if every sample ≤ 1 tip
}

max_copies_per_sample <- function(tr) {
  if (!inherits(tr, "phylo")) return(NA_integer_)
  ids <- vapply(tr$tip.label, sample_from_label, character(1))
  ids <- ids[ids %in% all_samples]
  if (length(ids) == 0) return(0L)
  max(tabulate(match(ids, unique(ids))))
}
n_samples_with_dups <- function(tr) {
  if (!inherits(tr, "phylo")) return(NA_integer_)
  ids <- vapply(tr$tip.label, sample_from_label, character(1))
  ids <- ids[ids %in% all_samples]
  if (length(ids) == 0) return(0L)
  sum(table(ids) > 1)
}

analyze_gene_tree <- function(tr, og_id) {
  if (!inherits(tr, "phylo")) {
    return(data.frame(
      og_id, n_tips = 0, n_nigoni = 0, n_briggsae = 0,
      mono_nigoni = NA, mono_briggsae = NA,
      n_mixed_clades = NA, mixed_summary = NA_character_,
      is_single_copy = NA,
      # optional diagnostics:
      max_copies_per_sample = NA_integer_,
      n_samples_with_dups = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  
  sp_map <- species_of(tr$tip.label)
  n_nig  <- sum(sp_map == "nigoni",   na.rm = TRUE)
  n_bri  <- sum(sp_map == "briggsae", na.rm = TRUE)
  
  mono_nig <- is_mono(tr, "nigoni")
  mono_bri <- is_mono(tr, "briggsae")
  
  mixes <- mixed_clades(tr)
  mixed_summary <- if (length(mixes)) {
    paste0(
      vapply(mixes, function(m) {
        paste0("node", m$node, ": ",
               paste(names(m$species_counts), as.integer(m$species_counts), sep = "=", collapse = ","),
               " | samples{", paste(m$samples, collapse = ","), "}",
               " | tips[", paste(m$tips, collapse = ","), "]")
      }, character(1)),
      collapse = " || "
    )
  } else NA_character_
  
  data.frame(
    og_id = og_id,
    n_tips = length(tr$tip.label),
    n_nigoni = n_nig,
    n_briggsae = n_bri,
    mono_nigoni = mono_nig,
    mono_briggsae = mono_bri,
    n_mixed_clades = length(mixes),
    mixed_summary = mixed_summary,
    is_single_copy = is_single_copy_tree(tr),
    # optional diagnostics:
    max_copies_per_sample = max_copies_per_sample(tr),
    n_samples_with_dups = n_samples_with_dups(tr),
    stringsAsFactors = FALSE
  )
}

res2 <- do.call(
  rbind,
  Map(analyze_gene_tree, tr = trees, og_id = names(trees))
)

res3 <- res2 %>% 
  dplyr::filter(is_single_copy==T) %>% 
  #dplyr::select(-mixed_summary) %>%
  dplyr::filter(mono_nigoni==F | mono_briggsae==F)



BL_THR      <- 1e-6   # ultra-short branches are "weak"
MIN_STRONG  <- 2      # >= this many strong internal edges -> structured
collapse_weak <- function(tr, bl_thr = BL_THR) {
  if (!inherits(tr, "phylo")) return(tr)
  te <- tr
  if (is.null(te$edge.length)) te$edge.length <- rep(1, nrow(te$edge))
  # mark edges shorter than threshold as weak
  weak <- te$edge.length < bl_thr
  te$edge.length[weak] <- 0
  di2multi(te, tol = 0)
}

# Count "strong" internal edges using ONLY branch length
count_strong_edges <- function(tr, bl_thr = BL_THR) {
  if (!inherits(tr, "phylo") || tr$Nnode == 0) return(0L)
  ne <- Ntip(tr)
  internal_idx <- which(tr$edge[,2] > ne)
  if (length(internal_idx) == 0 || is.null(tr$edge.length)) return(0L)
  sum(tr$edge.length[internal_idx] >= bl_thr)
}

# Monophyly after BL-only collapse (keeps your is_mono/species_of)
mono_after_collapse <- function(tr) {
  tc <- collapse_weak(tr)  # uses BL_THR only
  c(
    mono_nigoni_collapsed   = is_mono(tc, "nigoni"),
    mono_briggsae_collapsed = is_mono(tc, "briggsae")
  )
}


assess_og_bl <- function(og_id,
                         collapse_q = 0.05,   # collapse bottom 5% internal BLs
                         long_q = 0.25,       # "long" edges = top 75% internal BLs
                         min_long_edges = 2,  # need at least this many to call structured
                         min_internal_nodes = 3) {
  tr <- trees[[og_id]]
  if (!inherits(tr, "phylo")) {
    return(data.frame(
      og_id, mono_nigoni_collapsed = NA, mono_briggsae_collapsed = NA,
      n_long_edges = NA_integer_, n_internal_after = NA_integer_,
      polyphyly_class = "ambiguous", uninformative_flag = TRUE, stringsAsFactors = FALSE
    ))
  }

  tc   <- collapse_weak_bl(tr, collapse_q = collapse_q)
  mons <- mono_after_collapse_bl(tr, collapse_q = collapse_q)
  still_poly <- (isFALSE(mons[["mono_nigoni_collapsed"]]) || isFALSE(mons[["mono_briggsae_collapsed"]]))

  n_long <- count_long_edges_bl(tc, long_q = long_q)
  nin    <- if (inherits(tc, "phylo")) tc$Nnode else NA_integer_

  polyphyly_class <- if (!still_poly) {
    "star_like"                         # polyphyly vanished after collapse
  } else if (!is.na(n_long) && n_long >= min_long_edges) {
    "structured"                        # still poly + enough long internal edges
  } else if (!is.na(nin) && nin >= min_internal_nodes) {
    "structured"                        # still poly + enough resolution overall
  } else {
    "ambiguous"
  }

  data.frame(
    og_id = og_id,
    mono_nigoni_collapsed = mons[["mono_nigoni_collapsed"]],
    mono_briggsae_collapsed = mons[["mono_briggsae_collapsed"]],
    n_long_edges = n_long,
    n_internal_after = nin,
    polyphyly_class = polyphyly_class,
    uninformative_flag = polyphyly_class != "structured",
    stringsAsFactors = FALSE
  )
}

assess_og <- function(og_id) {
  tr <- trees[[og_id]]
  if (!inherits(tr, "phylo")) {
    return(data.frame(
      og_id, mono_nigoni_collapsed = NA, mono_briggsae_collapsed = NA,
      n_strong_edges = NA_integer_, polyphyly_class = "ambiguous",
      uninformative_flag = TRUE, stringsAsFactors = FALSE
    ))
  }
  mons <- mono_after_collapse(tr)
  tc   <- collapse_weak(tr)
  nse  <- count_strong_edges(tc)
  
  # classify: if polyphyly disappears (both not FALSE) -> star_like
  # else if still polyphyletic but has structure -> structured
  # else ambiguous
  still_poly <- (isFALSE(mons[["mono_nigoni_collapsed"]]) || isFALSE(mons[["mono_briggsae_collapsed"]]))
  polyphyly_class <- if (!still_poly) {
    "star_like"
  } else if (nse >= MIN_STRONG) {
    "structured"
  } else {
    "ambiguous"
  }
  
  # flag potentially uninformative (= not "structured")
  uninformative_flag <- polyphyly_class != "structured"
  
  data.frame(
    og_id = og_id,
    mono_nigoni_collapsed = mons[["mono_nigoni_collapsed"]],
    mono_briggsae_collapsed = mons[["mono_briggsae_collapsed"]],
    n_strong_edges = nse,
    polyphyly_class = polyphyly_class,
    uninformative_flag = uninformative_flag,
    stringsAsFactors = FALSE
  )
}
# Re-run assessment on your filtered set (res3) with BL-only logic
assessed <- do.call(rbind, lapply(res3$og_id, assess_og))

res3_flagged <- res3 %>%
  left_join(assessed, by = "og_id") %>%
  dplyr::select(-mixed_summary) %>%
  dplyr::filter(polyphyly_class=="structured")


extract_ref_gene <- function(tr) {
  if (!inherits(tr, "phylo")) return(NA_character_)
  hits <- grep("^QX1410_", tr$tip.label, value = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits[1]  # return the first (there should only be one)
}

# Apply to all trees in your list
ref_genes <- vapply(trees, extract_ref_gene, character(1))

res3_flagged$ref_gene <- ref_genes[res3_flagged$og_id]

gff <- readr::read_tsv("../processed_data/c_briggsae.QX1410_20250929.csq.gff",col_names = c("CHROM","source","type","start","end","score","strand","phase","attributes"))
tran <- gff %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("name","rest0"),sep=";Parent=") %>%
  tidyr::separate(rest0,into=c("rest1","loc"),sep=";locus=") %>%
  tidyr::separate(loc,into=c("locus","rest2"),sep=";") %>%
  tidyr::separate(rest1,into=c("parent","rest3"),sep=";biotype=") %>%
  dplyr::mutate(name=gsub("ID=transcript:","",name),parent=gsub("gene:","",parent)) %>%
  tidyr::separate(parent,into=c("parent_name","sequence_name"),sep=";") %>%
  dplyr::mutate(sequence_name=gsub("sequence_name=","",sequence_name)) %>%
  dplyr::select(CHROM,start,end,strand,name,sequence_name) 

res3_QX <- res3_flagged %>%
  dplyr::filter(!is.na(ref_gene)) %>%
  dplyr::mutate(ref_gene=gsub("QX1410_","",ref_gene)) %>%
  dplyr::left_join(tran,by=c("ref_gene"="name"))


hdrs <- readr::read_tsv("../processed_data/HDR_CB_allStrain_5kbclust_20250930.tsv") %>%
  dplyr::filter(STRAIN %in% all_str)


genes <- as.data.table(res3_QX)[, .(og_id, CHROM, start, end)]
hdrs_dt <- as.data.table(hdrs)[, .(CHROM, minStart, maxEnd, STRAIN)]

genes_pos <- genes[, .(og_id, CHROM, pos_start = start, pos_end = start)]
setkey(hdrs_dt, CHROM, minStart, maxEnd)
setkey(genes_pos, CHROM, pos_start, pos_end)

# overlap: gene point "within" HDR interval
ov <- foverlaps(genes_pos, hdrs_dt,
                by.x = c("CHROM","pos_start","pos_end"),
                type = "any", nomatch = 0L)

# summarise per gene
ov_sum <- ov[, .(
  in_hdr = TRUE,
  hdr_strains = paste(sort(unique(STRAIN)), collapse = ";"),
  n_hdr_hits = .N
), by = og_id]

res3_QX_annot <- res3_QX %>% 
  dplyr::left_join(as.data.frame(ov_sum),by="og_id") %>%
  dplyr::filter(!is.na(in_hdr))


hdr_tip_labels <- function(gt, hdr_strains_str) {
  if (is.null(hdr_strains_str) || is.na(hdr_strains_str) || nchar(hdr_strains_str) == 0) return(character(0))
  toks <- str_split(hdr_strains_str, ";")[[1]] |> trimws()
  toks <- toks[nzchar(toks)]
  if (!length(toks)) return(character(0))
  gt$tip.label[vapply(gt$tip.label, function(lbl) any(stringr::str_detect(lbl, fixed(toks))), logical(1))]
}

extract_tip_isotype_token <- function(lbl) {
  m <- stringr::str_match(lbl, "^([^_.-]+)")
  m[,2]
}

# 2) Robust mapper: prefer exact token match; else longest whole-word match across iso_vec
label_to_isotype_from_lineages <- function(lbl, iso_vec) {
  tok <- extract_tip_isotype_token(lbl)
  if (!is.na(tok) && tok %in% iso_vec) return(tok)
  
  # fallback: whole-word (non-alnum-bounded) matches; choose the longest to avoid NIC21 vs NIC2150 issues
  hits <- iso_vec[vapply(
    iso_vec,
    function(x) stringr::str_detect(lbl, regex(paste0("(?<![A-Za-z0-9])", x, "(?![A-Za-z0-9])"))),
    logical(1)
  )]
  if (length(hits)) return(hits[which.max(nchar(hits))])
  
  NA_character_
}



ogs_to_plot <- res3_QX_annot$og_id

x_max <- max(
  map_dbl(ogs_to_plot, ~ {
    tr <- trees[[.x]]
    if (is.null(tr)) NA_real_ else max(node.depth.edgelength(tr))
  }),
  na.rm = TRUE
)

make_tree_plot <- function(pick, trees, res3_QX_annot, label_to_species, lineages,
                           x_max = NULL, pad_frac = 0.25, lab_frac = 0.02) {
  gt <- trees[[pick]]
  if (is.null(gt)) return(NULL)
  this_row <- dplyr::filter(res3_QX_annot, og_id == pick)
  
  # species + lineage metadata per tip
  iso_vec <- lineages$isotype
  tip_isotype <- vapply(gt$tip.label, label_to_isotype_from_lineages, character(1), iso_vec = iso_vec)
  
  tip_meta <- tibble(
    label   = gt$tip.label,
    species = vapply(gt$tip.label, label_to_species, character(1)),
    isotype = tip_isotype
  ) %>%
    mutate(is_hdr  = label %in% hdr_tip_labels(gt, this_row$hdr_strains)) %>%
    left_join(lineages %>% dplyr::select(isotype, Lineage, lineage_color), by = "isotype") %>%
    mutate(
      Lineage       = if_else(is.na(Lineage) | !nzchar(Lineage), "Unknown", Lineage),
      lineage_color = if_else(is.na(lineage_color) | !nzchar(lineage_color), "grey60", lineage_color),
      is_hdr        = ifelse(species == "nigoni", "Unknown", is_hdr)
    )
  
  # palettes
  species_pal <- c(nigoni = "tomato", briggsae = "steelblue", unknown = "grey50")
  present_lineages <- unique(tip_meta$Lineage)
  palette_tbl <- lineages %>% dplyr::filter(Lineage %in% present_lineages) %>% dplyr::distinct(Lineage, lineage_color)
  lineage_pal <- setNames(rep("grey60", length(present_lineages)), present_lineages)
  if (nrow(palette_tbl)) lineage_pal[palette_tbl$Lineage] <- palette_tbl$lineage_color
  lineage_pal["Unknown"] <- "grey60"
  
  # --- axis scaling: use global x_max if provided; else per-tree span
  tree_span <- max(node.depth.edgelength(gt))
  span_ref  <- if (is.null(x_max)) tree_span else x_max
  right_pad <- span_ref * pad_frac
  lab_off   <- span_ref * lab_frac
  
  p <- ggtree(gt) %<+% tip_meta +
    geom_tiplab(aes(color = species), size = 3,
                offset = lab_off, align = TRUE, linesize = 0.2, linetype = "dotted") +
    scale_color_manual(values = c(nigoni = "tomato", briggsae = "steelblue", unknown = "grey50"),
                       guide = "none") +
    ggnewscale::new_scale("fill") +
    geom_tippoint(aes(subset = isTip, shape = is_hdr, fill = species),
                  color = "black", size = 3.2, stroke = 1) +
    scale_fill_manual(
      name = "Species",
      values = c(nigoni = "tomato", briggsae = "steelblue", unknown = "grey50"),
      breaks = c("briggsae", "nigoni", "unknown"),
      labels = c(expression(italic("C. briggsae")),
                 expression(italic("C. nigoni")),
                 "Unknown"),
      guide = guide_legend(override.aes = list(shape = 21, color = "black"))
    ) +
    scale_shape_manual(
      name = "HDR",
      breaks = c(FALSE, TRUE, "Unknown"),
      values = c(`FALSE` = 21, `TRUE` = 22, `Unknown` = 23),
      labels = c("No", "Yes", "Unknown")
    ) +
    theme_tree2() +
    ggtitle(
      bquote(
        .(this_row$og_id) * " - " * italic(.(this_row$sequence_name)) *
          " - " * .(this_row$CHROM) * ":" * .(this_row$start) * "-" * .(this_row$end)
      )
    ) +
    scale_x_continuous(limits = c(0, span_ref + right_pad),
                       expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(clip = "off")
  
  p
}

lineages <- readr::read_tsv("../processed_data/isotype_byLineage_GeoLocAdmCol_20250909.tsv")
# rebuild the list (now passing `lineages`)
plots_list <- setNames(
  lapply(
    res3_QX_annot$og_id,
    make_tree_plot,
    trees = trees,
    res3_QX_annot = res3_QX_annot,
    label_to_species = label_to_species,
    lineages = lineages,
    x_max = x_max,        # <- global max controls every plot's x-scale
    pad_frac = 0.25,      # consistent right padding for labels
    lab_frac = 0.02       # consistent label offset
  ),
  res3_QX_annot$og_id
)

# 
# test <- setNames(lapply(
#   "OG0000001",
#   make_tree_plot,
#   trees = trees,
#   res3_QX_annot = res3_QX_annot,
#   label_to_species = label_to_species,
#   lineages = lineages,
#   #x_max = NULL,        # <- global max controls every plot's x-scale
#   pad_frac = 0.25,      # consistent right padding for labels
#   lab_frac = 0.02       # consistent label offset
# ),
# "OG0000001"
# )

write.table(res3_QX_annot,"../tables/TableS11_introgression_trees.tsv",quote = F,col.names = F,row.names = F,sep = "\t")

label_has_any <- function(label, vec) {
  if (length(vec) == 0) return(FALSE)
  any(vapply(vec, function(p) grepl(p, label, fixed = TRUE), logical(1)))
}

# count how many of 'labels' contain ANY of the nigoni patterns
count_nigoni_labels <- function(labels, nigoni_vec) {
  sum(vapply(labels, label_has_any, logical(1), vec = nigoni_vec))
}

# for a tree 'gt' and a tip label 'tip_lbl', return the k nearest neighbor tip labels (excluding self)
k_nearest_neighbors <- function(gt, tip_lbl, k = 20) {
  D <- cophenetic(gt)                        # patristic distances
  if (!(tip_lbl %in% rownames(D))) return(character(0))
  drow <- D[tip_lbl, ]
  drow <- drow[names(drow) != tip_lbl]       # drop self
  ord  <- order(drow, decreasing = FALSE)
  nn   <- names(drow)[ord]
  head(nn, k)
}

# Given one row (og_id) compute the 3 outputs
summarize_og_neighbors <- function(og_id, trees, hdr_strains_str, nigoni_vec) {
  gt <- trees[[og_id]]
  if (is.null(gt)) {
    return(tibble(
      og_id = og_id,
      nigoni_neighbors = FALSE,
      which_target = NA_character_,
      max_nn_count = 0L
    ))
  }
  
  # Parse semicolon-separated hdr_strains tokens
  tokens <- hdr_strains_str %>%
    { if (is.null(.) || is.na(.) || !nzchar(.)) "" else . } %>%
    str_split(";", simplify = FALSE) %>%
    pluck(1) %>%
    trimws() %>%
    discard(~ !nzchar(.))
  
  if (length(tokens) == 0) {
    return(tibble(
      og_id = og_id,
      nigoni_neighbors = FALSE,
      which_target = NA_character_,
      max_nn_count = 0L
    ))
  }
  
  tip_labels <- gt$tip.label
  
  # For each token, find matching tip(s) by substring, then get k-NN and count nigoni neighbors
  per_token <- map(tokens, function(tok) {
    matched_tips <- tip_labels[stringr::str_detect(tip_labels, fixed(tok))]
    if (length(matched_tips) == 0) {
      return(tibble(token = tok, max_count = 0L))
    }
    # a token can match multiple tips; take the max count across those matched tips
    counts <- vapply(matched_tips, function(tip_lbl) {
      nn <- k_nearest_neighbors(gt, tip_lbl, k = 20)
      count_nigoni_labels(nn, nigoni_vec)
    }, integer(1))
    tibble(token = tok, max_count = max(counts))
  }) %>% list_rbind()
  
  any_hits     <- any(per_token$max_count > 0)
  which_tokens <- per_token %>% filter(max_count > 0) %>% pull(token)
  which_str    <- if (length(which_tokens)) paste(which_tokens, collapse = ";") else NA_character_
  max_count    <- if (nrow(per_token)) max(per_token$max_count) else 0L
  
  tibble(
    og_id = og_id,
    nigoni_neighbors = any_hits,
    which_target = which_str,
    max_nn_count = as.integer(max_count)
  )
}

# ---------- run over all og_id and merge back ----------

# 'nigoni' should be your character vector of substrings/IDs defining C. nigoni strains.
# e.g., nigoni <- c("NIC1593","QG3661", ...)  # <- you already have this
# 'trees' is your named list of phylo objects; names should include the og_id keys
# 'res3_QX_annot' is your data frame with og_id and hdr_strains

neighbor_summary <- map_dfr(
  res3_QX_annot$og_id,
  ~ summarize_og_neighbors(.x, trees, res3_QX_annot$hdr_strains[match(.x, res3_QX_annot$og_id)], nigoni)
)

res3_QX_annot_aug <- res3_QX_annot %>%
  left_join(neighbor_summary, by = "og_id") %>%
  dplyr::filter(!is.na(which_target))

plots_list_aug <- setNames(
  lapply(
    res3_QX_annot_aug$og_id,
    make_tree_plot,
    trees = trees,
    res3_QX_annot = res3_QX_annot_aug,
    label_to_species = label_to_species,
    lineages = lineages
  ),
  res3_QX_annot_aug$og_id
)


pimm_dir <- "../processed_data/intro_pimm"
  
ogs <- res3_QX_annot_aug$og_id
pimm_list <- setNames(vector("list", length(ogs)), ogs)

for (og in ogs) {
  file_path <- file.path(pimm_dir, paste0(og, ".pimm"))
  
  if (file.exists(file_path)) {
    pimm_list[[og]] <- read.csv(file_path, row.names = 1, check.names = FALSE)
  } else {
    warning(paste("File not found:", file_path))
    pimm_list[[og]] <- NULL
  }
}

make_heatmap <- function(og, pimm_list, trees) {
  
  
  gt <- trees[[og]]
  if (is.null(gt)) return(NULL)
  this_row <- dplyr::filter(res3_QX_annot, og_id == og)
  
  # species + lineage metadata per tip
  iso_vec <- lineages$isotype
  tip_isotype <- vapply(gt$tip.label, label_to_isotype_from_lineages, character(1), iso_vec = iso_vec)
  
  tip_meta <- tibble(
    label   = gt$tip.label,
    species = vapply(gt$tip.label, label_to_species, character(1)),
    isotype = tip_isotype
  ) %>%
    mutate(
      #species = ifelse(is.na(species), "unknown", species),
      is_hdr  = label %in% hdr_tip_labels(gt, this_row$hdr_strains)
    ) %>%
    left_join(lineages %>% select(isotype, Lineage, lineage_color), by = "isotype") %>%
    mutate(
      Lineage       = if_else(is.na(Lineage) | !nzchar(Lineage), "Unknown", Lineage),
      lineage_color = if_else(is.na(lineage_color) | !nzchar(lineage_color), "grey60", lineage_color)
    ) %>% dplyr::mutate(is_hdr=ifelse(species=="nigoni","NA",is_hdr))
  
  mat <- pimm_list[[og]]
  tr  <- trees[[og]]
  
  if (is.null(mat) || is.null(tr)) {
    warning(paste("Missing data for", og))
    return(NULL)
  }
  
  if (inherits(mat, "data.frame")) {
    mat <- as.matrix(mat)
  }
  
  if (!is.numeric(mat)) {
    stop("Matrix for ", og, " is not numeric after coercion.")
  }
  
  # ---- Clean sequence names ----
  clean_names <- function(x) {
    gsub("_longest_prot(_[Tt]ranscript)?", "", x)
  }
  
  # Clean matrix row/colnames
  rownames(mat) <- clean_names(rownames(mat))
  colnames(mat) <- clean_names(colnames(mat))
  
  asym <- max(abs(mat - t(mat)), na.rm = TRUE)
  if (is.finite(asym) && asym > 1e-8) {
    message("Matrix ", og, " is asymmetric (max |A - A^T| = ", signif(asym, 3), ").")
  }
  
  if (any(duplicated(rownames(mat))) || any(duplicated(colnames(mat)))) {
  stop("Duplicate row/column names after cleaning in ", og,
       ". Consider disambiguating names before building the matrix.")
}
  
  # Clean tree tip labels
  tr$tip.label <- clean_names(tr$tip.label)
  
  # ---- Match taxa ----
  common <- intersect(rownames(mat), tr$tip.label)
  if (length(common) < 3) {
    warning(paste("Too few common taxa in", og))
    return(NULL)
  }
  
  mat <- mat[common, common, drop = FALSE]
  tr  <- keep.tip(tr, common)
  
  # Order matrix by tree tip order
  tip_order <- tr$tip.label
  mat <- mat[tip_order, tip_order]
  
  # set clustered dendrogram 
  d <- as.dist(1 - mat)         # symmetric if mat is symmetric
  hc <- hclust(d, method = "average")
  dend <- as.dendrogram(hc)
  
  tags <- sapply(strsplit(rownames(mat), "_"), `[`, 1)
  row_colors <- ifelse(tags %in% nigoni, "tomato",
                       ifelse(tags %in% briggsae, "steelblue", "gray"))
  # ---- Define color scale ----
  col_fun <- colorRamp2(
    c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    c("white", "white", "lightblue", "blue", "yellow", "orange", "red")
  )
  
  
  tip_meta <- tip_meta %>%
    filter(label %in% tip_order) %>%
    slice(match(tip_order, label))
  
  # Convert is_hdr (TRUE/FALSE/NA) to a factor for consistent coloring
  tip_meta$is_hdr <- factor(
    tip_meta$is_hdr,
    levels = c("TRUE", "FALSE", "Unknown"),
    labels = c("HDR", "Non-HDR", "Unknown")
  )
  
  # Define colors for annotation
  hdr_colors <- c(
    "HDR" = "red",
    "Non-HDR" = "black",
    "NA" = "grey80"
  )
  
  species_colors <- c(
    briggsae = "steelblue",
    nigoni   = "tomato"
  )
  
  # Create the bottom column annotation
  bottom_anno <- HeatmapAnnotation(
    HDR     = tip_meta$is_hdr,     # factor with levels/labels you set earlier
    species = factor(tip_meta$species, levels = c("briggsae","nigoni")),
    col = list(
      HDR     = hdr_colors,
      species = species_colors
    ),
    annotation_name_side = "left",
    annotation_legend_param = list(
      HDR = list(
        title  = "HDR",
        at     = c("HDR", "Non-HDR", "NA"),
        labels = c("Yes", "No", "NA")
      ),
      species = list(
        title  = "Species",
        at     = c("briggsae","nigoni"),
        # italicized genus/species in the legend:
        labels = c(
          expression(italic("C. briggsae")),
          expression(italic("C. nigoni"))
        )
      )
    )
  )
  
  # ---- Create Heatmap ----
  hm <- grid::grid.grabExpr(draw(Heatmap(
      mat,
      name = "Identity",
      col = col_fun,
      cluster_rows = dend,
      cluster_columns = dend,
      # clustering_distance_rows = function(m) as.dist(1 - m),
      # clustering_distance_columns = function(m) as.dist(1 - m),
      # clustering_method_rows = "average",
      # clustering_method_columns = "average",
      show_row_names = TRUE,
      show_column_names = FALSE,
      row_names_gp = grid::gpar(col = row_colors, fontsize = 9),
      heatmap_legend_param = list(
        title = "A.A.\n%ID",
        at = seq(0.5, 1.0, by = 0.1),
        labels = seq(50, 100, by = 10)  # <- multiplied by 100
      ),
      bottom_annotation = bottom_anno
      ),
    merge_legend = TRUE,
    )
  )
  return(hm)
}
ogs <- names(pimm_list)

heatmap_list <- map(ogs, ~ make_heatmap(.x, pimm_list, trees)) %>%
  set_names(ogs)

ht_to_ggplot <- function(ht) {
  #g <- grid::grid.grabExpr(ComplexHeatmap::draw(ht), wrap = TRUE)
  ggplotify::as.ggplot(ht)
}

combined_plots <- list()

for (og in names(heatmap_list)) {
  if (og %in% names(plots_list_aug)) {
    ht_gg <- ht_to_ggplot(heatmap_list[[og]])  # ggplotified heatmap
    gp    <- plots_list_aug[[og]]              # your ggplot
    
    combined_plots[[og]] <- cowplot::plot_grid(
      gp, ht_gg, 
      ncol = 1  # use nrow = 2 for vertical stacking
    )
  }
}


out_dir <- "../figures/FileS2_introgression_figures/"   
#dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# lookup table for filename parts
fn_lookup <- res3_QX_annot_aug %>%
  select(og_id, max_nn_count, sequence_name)

# sanitize to safe filenames
sanitize <- function(x) {
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")  # replace unsafe chars with _
  x <- str_squish(x)
  x
}

# save all plots
iwalk(combined_plots, function(p, og) {
  if (is.null(p) || !inherits(p, "ggplot")) return(invisible(NULL))
  row <- fn_lookup %>% filter(og_id == og) %>% slice(1)
  
  maxNN <- if (nrow(row)) row$max_nn_count else NA_integer_
  seqnm <- if (nrow(row)) row$sequence_name else "unknown"
  
  # prefix = numeric maxNN (or "NA" if missing)
  prefix <- if (is.na(maxNN)) "NA" else as.character(maxNN)   # or formatC(maxNN, digits=0)
  
  # 12maxNN_OGXXXX_sequence_name.png
  fname <- paste0(
    prefix, "maxNN_", sanitize(og), "_", sanitize(seqnm), ".png"
  )
  
  # optional: keep filenames reasonably short
  if (nchar(fname) > 180) {
    base <- tools::file_path_sans_ext(fname)
    ext  <- tools::file_ext(fname)
    base <- substr(base, 1, 170)
    fname <- paste0(base, ".", ext)
  }
  
  ggsave(
    filename = file.path(out_dir, fname),
    plot     = p,
    width    = 7, height = 9, units = "in", dpi = 600, bg = "white"
  )
})

draw_colored_tree_ggtree <- function(
    gt,
    label_to_species,
    palette = c(nigoni = "tomato", briggsae = "steelblue", unknown = "grey50"),
    title  = "Concordant",
    label_size = 3,
    offset = 0.02,
    align = TRUE
) {
  stopifnot(!is.null(gt), !is.null(gt$tip.label))
  
  # --- keep original labels for species mapping
  orig_labels <- gt$tip.label
  
  # map labels -> species (using ORIGINAL labels)
  sp <- vapply(orig_labels, label_to_species, character(1))
  sp[is.na(sp) | !nzchar(sp)] <- "unknown"
  
  # --- create a DISPLAY label (cleaned AFTER species is set)
  label_disp <- sub("^QX1410_", "", orig_labels)       # remove leading prefix if present
  label_disp <- sub("\\.t[0-9]+$", "", label_disp)     # remove trailing .t#
  
  # ensure unknown color exists
  if (!"unknown" %in% names(palette)) palette <- c(palette, unknown = "grey50")
  
  # metadata joined by 'label' (must match gt$tip.label)
  tip_meta <- tibble::tibble(label = orig_labels, species = sp, label_disp = label_disp)
  
  # dynamic offset
  span <- max(ape::node.depth.edgelength(gt))
  off  <- if (offset < 1) span * offset else offset
  
  p <- ggtree::ggtree(gt) %<+% tip_meta +
    ggtree::geom_tiplab(
      ggplot2::aes(color = species, label = label_disp),
      size = label_size,
      offset = off,
      align = align,
      linesize = if (align) 0.2 else NA,
      linetype = if (align) "dotted" else "solid"
    ) +
    ggplot2::scale_color_manual(
      values = palette,
      breaks = c("briggsae", "nigoni", "unknown"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, span + off * 3)) +
    ggtree::theme_tree2() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = grid::unit(c(0, 5, 0, 5), "mm")) +
    ggplot2::ggtitle(title)
  
  return(p)
}
pick <- "OG0012935"
gt   <- trees[[pick]]

pick2 <- "OG0018191"
gt2 <- trees[[pick2]]

cc_tree <- draw_colored_tree_ggtree(
  gt,
  label_to_species = label_to_species,
  title = "",
  label_size = 2,
  offset = 0.2, 
  align = TRUE
)

dc_tree <- draw_colored_tree_ggtree(
  gt2,
  label_to_species = label_to_species,
  title = "",
  label_size = 2,
  offset = 0.2,   
  align = TRUE
)


consensus_plt <- draw_colored_tree_ggtree(
  consensus,
  label_to_species = label_to_species,
  title = "",
  label_size = 2,
  offset = 0.05,  
  align = TRUE
)

S18 <- cowplot::plot_grid(consensus_plt,cc_tree,dc_tree,nrow=1,labels = c("a","b","c"))

ggsave(S18,filename = "../figures/Figure_S18_treesamples.png",width    = 7, height = 4, units = "in", dpi = 600, bg = "white")
