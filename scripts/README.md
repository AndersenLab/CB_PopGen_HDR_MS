# this readme is currently under development

# Scripts
This directory contains the scripts used to carry out analysis and generate figures and tables for the manuscript.

## Directory Structure

scripts/
- geography_trees_PCA
- genetic_similarity_and_admixture/
- diversity_and_divergence/
- gene_enrichment/
- HDRs/

## Abreviations
E.D.F. = Extended Data Figure
S.F. = Supplementary Figure
S.D. = Supplementary Data

## geography_trees_PCA

- `Calculate_Cb_generate_link_and_sample_list.sh`
script description
Generates?

- `Calculate_Cb_pyh_to_tree.sh`  
script description
Generates?

- `Calculate_Cb_vcf_to_pyh.sh`  
script description
Generates?

- `Calculate_Cb_pca_by_chrom.sh`  
script description
Generates?

- `Cb_PCA_by_chrom_LD_0.9.R`  
script description
Generates Figure S5?

- `nf_Cb_pca_by_chrom.sh`
script description
Generates?

- `Cb_iterations_5_10_15_20_25_30_35_40.sh`  
script description
Generates?

- `Cb_par_10_iterations.par`  
script description
Generates?

- `Cb_par_15_iterations.par`
script description
Generates?

- `Cb_par_20_iterations.par`  
script description
Generates?

- `Cb_par_25_iterations.par`
script description
Generates?

- `Cb_par_30_iterations.par`
script description
Generates?

- `Cb_par_35_iterations.par`
script description
Generates?

- `Cb_par_40_iterations.par`
script description
Generates?

- `Cb_par_5_iterations.par`
script description
Generates?

- `Figure_S6_Cb_all_isotypes_PCA_LD_0.9_5_10_15_iterations.R`
script description
Generates Figure S5 ?

- `Cb_all_isotypes_PCA_LD_0.9.R`
script description
Generates Table S3 ? 

- `Cb_assemble_figure_1.R` 
script description
Generates Figure 1 ?

`Cb_Geo_locations_strains.R` 
script description
Generates?

`nf_Cb_pruned_VCF_and_PCA.sh`  
script description
Generates?

`Plot_tree_with_previous_lineages.R`
script description
Generates Figure S4 ?

- `Cb_Geo_locations_isotype.R`  
script description
Generates?

- `Figure_S1.R`  
script description
Generates Figure S1?

## genetic_similarity_and_admixture

- `Define_cosmopolitan_isotypes.R`  
script description
Generates?

- `Figure_S2_Plot_cosmopolitan_17_heatmap.R`
script description
Generates Figure S2?

- `Cb_concordance_histogram_hard_filtered.R`
script description
Generates Figure S23

- `Figure_S28_Cb_concordance_histogram_hard_filtered.R`
script description
I think this script is obsolete (replaced by `Cb_concordance_histogram_hard_filtered.R`?)

- `Figure_S3_Plot_cosmopolitan_concordance.R`
script description
Generates Figure S3

- `Calculate_Cb_admixture_post_processing.sh`  
script description
Generates?

- `Calculate_Cb_run_admixture_arrays.sh`  
script description
Generates Figure S9?

- `Plot_Cb_admixture_by_lineage.R`  
script description
Generates Figure S10 and S11

- `Plot_Cb_admixture_CV.R`
script description
I think this script is obsolete (replaced by `genetic_similarity_admixture.R`?)

- `Table_S2_concordance_strains.R`
script description
I think this script is obsolete (Table S2 no longer exists. Pairwise concordance is present in `processed_data`)

- `genetic_similarity_admixture.R`
script description
Generates Figure 2, S7, S8

## diversity_and_divergence

- `Cb_vcf_to_zarr.sh`
script description
Generates?

- `Cb_zarr_to_pi_theta_d_10kb.sh`
script description
Generates?

- `nucleotide_diversity_and_divergence.R`
script description
Generates Figure 3, 5, S12, S21, S22

- `Plot_pi_theta_d_lineage_three_figures.R`
script description
Generates Figures S19, S20

- `Table_Cb_geo_p_theta.R`
script description
Generates Table S5

- `Plot_Cb_pi_theta_d_geo_lineage_Autosomes_and_X.R`
script description
Generates Table S6

- `Table_S12_Plot_Cb_pi_theta_d_lineage_Autosomes_and_X.R`
script description
Generates Table S12

## gene_enrichment

- `InterProScan.sh`  
script description

- `IPR_GO_briggsae.R`
script description

Generates Figure 4, Table S9

## HDRs
- `call_HDRs.R`  
script description
Generates Figures S13, S25, S26, S27, S28, S29, S30, S31, and Table S7

- `characterize_HDRs.R`
script description
Generates Figures S14, S15, and Table S8

- `visualize_gene_content.R`
Given a reference region and a list of strains, generate gene content visualizations with orthologous relationships across wild strain genomes using genome alignments.
Generates sub-plots (as R objects) for merge_gene_content_plots.R
Genrates Figure S17

- `merge_gene_content_plots.R`
Given the paths to R objects from `visualize_gene_content.R`, generate Figure S16
  
