# this readme is currently under development

# Scripts
This directory contains the scripts used to carry out analysis and generate figures and tables for the manuscript.

## Directory Structure

scripts/
- process_dataset
- geography_trees_PCA
- genetic_similarity_and_admixture/
- diversity_and_divergence/
- gene_enrichment/
- HDRs/
- introgression_analysis/

## Abreviations
EDF = Extended Data Figure
SF = Supplementary Figure
SD = Supplementary Data

## process_dataset

script description

These scripts generate intermediate datasets reused across analyses and should be run before downstream scripts.

## geography_trees_PCA

- `Calculate_Cb_generate_link_and_sample_list.sh`

script description

Generates isotype list

- `Calculate_Cb_pyh_to_tree.sh` 

script description

Generates ML tree

- `Calculate_Cb_vcf_to_pyh.sh` 

script description

Generates PHYLIP format file for generating trees

- `Calculate_Cb_pca_by_chrom.sh` 

script description

Generates vcf files for each chromosome

- `Cb_PCA_by_chrom_LD_0.9.R`  

script description

Generates Supplementary Figure 2

- `nf_Cb_pca_by_chrom.sh`

script description

Generates PCA eigenvalues and PC scores for each chromosome

- `Cb_iterations_5_10_15_20_25_30_35_40.sh`  

script description

Generates PCA eigenvalues and PC scores after iterative outlier removal

- `SF3_Cb_all_isotypes_PCA_LD_0.9_5_10_15_iterations.R`

script description

Generates Supplementary Figure 3

- `Cb_all_isotypes_PCA_LD_0.9.R`

script description

Generates Supplementary Data 2

- `Cb_assemble_figure_1.R` 

script description

Generates Figure 1 

`Cb_Geo_locations_strains.R` 

script description

Generates panel a of Figure 1

`nf_Cb_pruned_VCF_and_PCA.sh`  

script description

Generates LD-pruned vcf file, PCA eigenvalues, and PC scores

`Plot_tree_with_previous_lineages.R`

script description

Generates Extended Data Figure 3 

- `Cb_Geo_locations_isotype.R`  

script description

Generates panel a of Extended Data Figure 1

- `make_EDF1.R`  

script description

Generates Extended Data Figure 1

## genetic_similarity_and_admixture

- `Define_cosmopolitan_isotypes.R`  

script description

Generates panel a of Supplementary Figure 1

- `SF1_Plot_cosmopolitan_17_heatmap.R`

script description

Generates Supplementary Figure 1

- `Cb_concordance_histogram_hard_filtered.R`

script description

Generates Supplementary Figure 14

- `EDF2_Plot_cosmopolitan_concordance.R`

script description

Generates Extended Data Figure 2

- `Calculate_Cb_admixture_post_processing.sh`  

script description

Generates summary table of ADMIXTURE CV errors

- `Calculate_Cb_run_admixture_arrays.sh`  

script description

Runs ADMIXTURE and generates Q and P matrices with CV results

- `Plot_Cb_admixture_by_lineage.R`  

script description

Generates S.F. 6 and S.F. 7

- `Plot_Cb_admixture_CV.R`

script description

Generates S.F. 5

- `genetic_similarity_admixture.R`

script description

Generates Figure 2

Generates Extended Data Figure 4

Generates Supplementary Figure 4

## diversity_and_divergence

- `Cb_vcf_to_zarr.sh`

script description

Converts the VCF files to Zarr files for downstream pi theta d analyses

- `Cb_zarr_to_pi_theta_d_10kb.sh`

script description

Calculates pi theta d from the Zarr files

- `pi_theta_d_python`

script description

Core scripts for calculating pi theta d

- `nucleotide_diversity_and_divergence.R`

script description

Generates Figures 3 and 5

Generates Supplementary Figure 8

Generates Extended Data Figures 8 and 9

- `Plot_pi_theta_d_lineage_three_figures.R`

script description

Generates Supplementary Figure 12

Generates Supplementary Figure 13 

- `Table_Cb_geo_p_theta.R`

script description

Generates Supplementary Data 4

- `Plot_Cb_pi_theta_d_geo_lineage_Autosomes_and_X.R`

script description

Generates Supplementary Data 5

- `SD12_Plot_Cb_pi_theta_d_lineage_Autosomes_and_X.R`

script description

Generates Supplementary Data 12

## gene_enrichment

- `InterProScan.sh`  

script description

Generates ? 

- `IPR_GO_briggsae.R`

script description

Generates Figure 4

Generates Supplementary Data 8

## HDRs
- `call_HDRs.R`  

Given variant counts and coverage statistics derived from long-read data  across reference genomic bins for select Tropical wild strains, optimize HDR calling parameters using short-read data for the Tropical relatedness group. Once optimal parameters are chosen, call HDRs across every relatedness group. The genomic positions of HDRs called in non-Tropical relatedness groups are transformed into the Tropical reference genome positions using genome-to-genome alignments.

Generates Supplementary Figures 9, 15, 16, 17, 18, 19, 20, 21

Generates Supplementary Data 6

- `characterize_HDRs.R`

Given variant counts and coverage statistics derived from short-read data across reference genomic bins for all wild strains, and a set of hyper-divergent regions called relative to QX1410 or a relatedness group reference genome, estimate the proportion of variants within HDRs and the proportion of the genome covered by HDRs.

Generates Supplementary Figure 10 

Generates Extended Data Figure 5 

Generates Supplementary Data 7

- `visualize_gene_content.R`

Given a reference region and a list of strains, generate gene content visualizations with orthologous relationships across wild strain genomes using genome alignments.

Generates sub-plots (as R objects) for `merge_gene_content_plots.R`.

Generates Supplementary Figure 11.

- `merge_gene_content_plots.R`

Merges and plots the R objects from `visualize_gene_content.R`.

Generates Extended Data Figure 6
  
## introgression_analysis

- `classify_tree_topology.R`

Given concatenated tree files of orthologous genes, identify single-copy ortholog trees with at least 1 branch lenth >1e6 where a reference gene is present and it is located within an HDR.

Generates Extended Data Figure 7

Generates  Supplementary Data 11

Visualizes tree figures in `~/figures/introgression_figures`