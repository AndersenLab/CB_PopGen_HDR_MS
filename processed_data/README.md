# Processed Data

This directory contains the processed datasets used for downstream population genomic analyses. Each subdirectory corresponds to a major analysis component. The files are derived from upstream pipelines (see `scripts/` and the main repository README for details).

## Directory Structure

processed_data/
- Geo_info/
- diversity_and_divergence/
- gene_diversity/
- genetic_similarity_and_admixutre/
- HDRs/
- gene_enrichment/

## Geo_info/

Processed files summarizing geographic info of the samples.

### Contains

- `Cb_cosmopolitan_isotype.txt`  
  List of C. briggsae isotypes classified as cosmopolitan.

- `Cb_geo_distance_within_same_isotype.csv`  
  Pairwise geographic distances among sampling locations within the same isotype.

- `Cb_indep_isotype_info_geo.csv`  
  Geographic information for each independent isotype.

- `Cb_indep_strain_info_geo.csv`  
  Geographic information for each independent strain with cosmopolitan strains labeled.

- `Cb_isotype_geo_freq.csv`  
  Frequency of each isotype across geographic regions, with cosmopolitan strains treated as a group.

- `indep_info_geo_for_each_strain.csv`  
  Consolidated geographic metadata for independent strains, providing specific location data for every strain, even though some of them are cosmopolitan.

- `strain_geo_freq.csv`  
  Frequency of each strain across geographic regions, providing summary data for every strain, even though some of them are cosmopolitan.

- `non_admixed_isotype_replicate_*.txt`  
  Each file corresponds to a replicate of the non-admixed isotype selection procedure.

## diversity_and_divergence/

Processed files summarizing genome-wide diversity and divergence metrics.

### Contains

- `chromosome_domain_Cbriggsae.csv`  
  Chromosomal domain breakpoints.

- `Dxy_Tropical/`  
  Pairwise absolute divergence (Dxy) between Tropical and other relatedness groups.

- `pi_theta_d/`  
  Nucleotide diversity metrics per genomic window.

- `pi_theta_d_geo/`  
  Nucleotide diversity metrics stratified by geographic grouping.

- `pi_theta_d_rg/`  
  Nucleotide diversity metrics stratified by relatedness group.

## gene_diversity/

Processed files used in gene-level visualization of HDRs.

### Contains

- `CBCN_L1L2_master.tsv`  [Large file. Available publicly on AWS] <br>
  L1 (gene) and L2 (mRNA) features from BRAKER gene predictions across all strains where a long-read genome was assembled.

- `CBCN_nucmer_db_20250603.tsv`  [Large file. Available publicly on AWS] <br>
  Whole-genome alignment (NUCmer) results.

- `CBCN_orthogroups.tsv`  [Large file. Available publicly on AWS] <br>
  Orthogroup assignments (Also includes _C. nigoni_ lines, but only _C. briggsae_ results are analyzed).

- `c_briggsae.QX1410_20250929.csq.gff`  [Large file. Available publicly on AWS] <br>
  Reference gene models.

- `HDR_I_*`, `HDR_II_*`, `HDR_V_*` (`.Rds`)  
  R objects containing gene-level diveristy plots for select hyper-divergent regions.

## genetic_similarity_and_admixutre/

Population structure, admixture, and genetic similarity data.

### Contains

- `best_k_long_admix_pops.csv`  
  Final selected ADMIXTURE results.

- `concat_Qfiles_K*.tsv`  
  Q-matrices for each value of K tested.

- `cv_matrix.tsv`  
  Cross-validation error matrix.

- `isotype_groups.tsv`  
  Isotype assignments.

- `non_admixed_isotypes.txt`  
  List of non-admixed strains.

- `K22_*` files (`.csv`, `.tsv`)  
  Best-K replicate results and processed ancestry tables.

- `gtcheck.tsv`
  Pairwise isotype genetic similarity estimates.

- `strain_gtcheck.tsv`  [Large file. Available publicly on AWS] <br>
  Pairwise strain genetic similarity estimates.

## HDRs/

Hyper-divergent region (HDR) calls associated and intermediate files.

### Contains

- `HDR_CB_allStrain_5kbclust_20250930.tsv`  
  HDR windows and across all strains relative to the Tropical (QX1410) reference genome.

- `HDR_CB_otherRG_UNT_5kbclust_20250930.tsv`  
  HDR windows for non-Tropical groups relative to their respective group's reference genome (UNT = untransformed).

- `Tropical.thresh_cov.tsv`  [Large file. Available publicly on AWS] <br>
  Coverage statistics for Tropical samples relative to the Tropical (QX1410) reference genome.

- `Tropical.variant_counts.tsv`  [Large file. Available publicly on AWS] <br>
  Variant counts computed for Tropical samples relative to the Tropical (QX1410) reference genome.

- `Other_RG.thresh_cov.tsv`  [Large file. Available publicly on AWS] <br>
  Coverage statistics for non-Tropical samples relative to their respective relatedness group reference genome.

- `Other_RG.variant_counts.tsv`  [Large file. Available publicly on AWS] <br>
  Variant counts for non-Tropical samples relative to their respective relatedness group reference genome.

- `CB_all_strain_vc.tsv`  [Large file. Available publicly on AWS] <br>
  Variant counts computed for all samples relative to the Tropical (QX1410) reference genome.

- `phy_file_LD_0.9.phy.contree`  
  Consensus phylogenetic tree of all isotypes.

- `QX1410_genomic_windows.1kb.bed`  
  BED file of Tropical (QX1410) reference 1 kb genomic windows.

- `relatednessGroup_nucmer_alignments.tsv`  
  Whole-genome alignment results between each relatedness group reference genomes and QX1410 (for mapping HDRs between relatedness groups).

- `Tropical_hifi_coords.tsv`  
  Whole-genome alignment results between each Tropical long-read genomes relative to QX1410 (for HDR parameter optimization).

- `Tropical_samples.txt`  
  List of tropical strains.

## gene_enrichment/

Data used for functional gene enrichment analysis to create Figure 4.

### Contains

- `c_briggsae.QX1410_20250929.csq.longest.nomito.prot.fa`
  QX1410 proteome created from the longest isoform transcripts and MtDNA removed (for InterProScan annotions)

- `c_elegans.N2.WS283.PConly.longestIsoform.gff3`
  C. elegans N2 GFF file for pulling N2 gene names from identified orthologs of QX1410 genes.

- `chromosome_domain_Cbriggsae.csv`
  QX1410 chromosome domains for classifying chromosomal arms.

- `QX1410_AF16_N2_orthogroups.tsv`
  OrthoFinder-identified orthologs among QX1410, AF16, and N2.

- `QX_IPR_allApps_20251006.tsv`
  Raw output from InterProScan.
