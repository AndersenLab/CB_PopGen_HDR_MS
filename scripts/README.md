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

These scripts generate intermediate datasets reused across analyses and should be run before downstream scripts.

## geography_trees_PCA

- `Calculate_Cb_generate_link_and_sample_list.sh`

Creates the isotype sample list used for downstream analyses.

- `Calculate_Cb_pyh_to_tree.sh` 

Generates maximum likelihood trees from PHYLIP-formatted alignments.

- `Calculate_Cb_vcf_to_pyh.sh` 

Generates PHYLIP format file for generating trees from VCF file

- `Calculate_Cb_pca_by_chrom.sh` 

Generates per-chromosome VCF files for plotting chromosome-specific PCA results.

Generates vcf files for each chromosome

- `Cb_PCA_by_chrom_LD_0.9.R`  

Plots chromosome-specific PCA results using LD-pruned variants.

Generates Supplementary Figure 2

- `nf_Cb_pca_by_chrom.sh`

Runs PCA on LD-pruned, per-chromosome VCF files and generates eigenvalues and PC scores.

Generates PCA eigenvalues and PC scores for each chromosome

- `Cb_iterations_5_10_15_20_25_30_35_40.sh`  

Performs iterative PCA analyses after sequential outlier removal and outputs eigenvalues and PC scores for each iteration.

Generates PCA eigenvalues and PC scores after iterative outlier removal

- `SF3_Cb_all_isotypes_PCA_LD_0.9_5_10_15_iterations.R`

Plots iterative PCA results across successive outlier-removal thresholds.

Generates Supplementary Figure 3

- `Cb_all_isotypes_PCA_LD_0.9.R`

Plots PCA results for all isotypes using LD-pruned variants.

Generates Supplementary Data 2

- `Cb_assemble_figure_1.R` 

Assembles panels and generates Figure 1.

Generates Figure 1 

`Cb_Geo_locations_strains.R` 

Plots geographic locations of strains for Figure 1

`nf_Cb_pruned_VCF_and_PCA.sh`  

Generates LD-pruned VCF files and performs PCA to obtain eigenvalues and PC scores.

Generates LD-pruned vcf file, PCA eigenvalues, and PC scores

`Plot_tree_with_previous_lineages.R`

Plots maximum likelihood tree annotated with previously defined lineages.

Generates Extended Data Figure 3 

- `Cb_Geo_locations_isotype.R`  

Plots geographic locations of isotype reference strains for Extended Data Figure 1.

- `make_EDF1.R`  

Assembles panels and generates Extended Data Figure 1.

Generates Extended Data Figure 1

## genetic_similarity_and_admixture

- `Define_cosmopolitan_isotypes.R`  

Defines cosmopolitan isotypes based on geographic distribution and generates summary classifications.

Generates panel a of Supplementary Figure 1

- `SF1_Plot_cosmopolitan_17_heatmap.R`

Plots genetic similarity heatmap of the 17 cosmopolitan isotypes.

Generates Supplementary Figure 1

- `Cb_concordance_histogram_hard_filtered.R`

Plots histogram of genetic similarity across hard-filtered variants.

Generates Supplementary Figure 14

- `EDF2_Plot_cosmopolitan_concordance.R`

Plots genetic similarity among cosmopolitan isotypes.

Generates Extended Data Figure 2

- `Calculate_Cb_admixture_post_processing.sh`  

Summarizes ADMIXTURE cross-validation errors across K values.

Generates summary table of ADMIXTURE CV errors

- `Calculate_Cb_run_admixture_arrays.sh`  

Runs ADMIXTURE across K values and generates Q and P matrices with cross-validation results.

- `Plot_Cb_admixture_by_lineage.R`  

Plots ADMIXTURE results by relatedness groups.

Generates Supplementary Figures 6 and 7

- `Plot_Cb_admixture_CV.R`

Plots ADMIXTURE cross-validation errors across K values.

Generates Supplementary Figure 5

- `genetic_similarity_admixture.R`

Using cross-validation ADMIXTURE results, select the optimal K parameter that minimizes shared subpopulation assignemnts between relatedness groups. Using pairwise genetic similarity estimates, visualize genetic similarity heatmap with strains ordered by hierarchical clustering with average linkage. Also generates heatmap of mean pairwise genetic similarity between relatedness groups.

Generates Figure 2

Generates Extended Data Figure 4

Generates Supplementary Figure 4

## diversity_and_divergence

- `Cb_vcf_to_zarr.sh`

Converts the VCF files to Zarr files for downstream pi theta d analyses

Generate Zarr files for calculating pi theta d.

- `Cb_zarr_to_pi_theta_d_10kb.sh`

Calculates pi, theta, and D in 10 kb windows from Zarr files.

Generate pi, theta, D summary table

- `pi_theta_d_python`

Core scripts for calculating pi theta d

- `nucleotide_diversity_and_divergence.R`

Visualizes bin-wise, genome-wide nucleotide diversity estimates (pi, thetaW) and physical position of hyper-divergent regions across all isotype reference strains. Classifies bin-wise Tajima's D and Dxy (for every pairwise relatedness group comparison) estimates based on overlaps with hyper-divergent regions. Visualizes bin-wise and mean Tajima's D and Dxy between hyper-divergent and non-hyper-divergent bins.

Generates Figures 3 and 5

Generates Supplementary Figure 8

Generates Extended Data Figures 8 and 9

- `Plot_pi_theta_d_lineage_three_figures.R`

Plots windowed pi and theta across chromosomes for all relatedness groups (facet: lineage and chromosome), with genome-domain backgrounds.

Generates Supplementary Figure 12

Generates Supplementary Figure 13 

- `Table_Cb_geo_p_theta.R`

Summarizes relatedness group-specific pi and theta.

Generates Supplementary Data 4

- `Plot_Cb_pi_theta_d_geo_lineage_Autosomes_and_X.R`

Summarizes mean windowed pi and theta for autosomes vs ChrX, split by chromosomal arms vs centers, across geographic/relatedness groups subsets.

Generates Supplementary Data 5

- `SD12_Plot_Cb_pi_theta_d_lineage_Autosomes_and_X.R`

Summarizes mean windowed pi and theta for autosomes and ChrX, split by chromosomal arms and centers, across relatedness groups.

Generates Supplementary Data 12

## gene_enrichment

- `InterProScan.sh`  

Script to annotate QX1410 protein sequences using InterProScan.

Generates InterProScan annotaitons table for `IPR_GO_briggsae.R`

- `IPR_GO_briggsae.R`

Given the QX1410 reference genome and InterProsScan annotations of the QX1410 proteome, this script classifies genes that are found in the chromosomal arm domains, genes found in hyper-divergent regions in chromosomal arm domains, and performs a one-sided hyper-geometric test to identify statistically-enriched InterProScan functional protein domains in genes found in hyper-divergent regions in chromosomal arm domains in relation to genes found outside of hyper-divergent regions in chromosomal arm domains.

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

Given concatenated tree files of orthologous genes, identify single-copy ortholog trees with at least 1 branch lenth >1e6 where a reference gene is present and it is located within an HDR. Visualize examples of concordant and discordant tree topologies relative to a consensus tree. Visualizes tree topologies that serve as examples of potential introgression, accompanied by amino-acid identity matrices.

Generates Extended Data Figure 7

Generates  Supplementary Data 11

Visualizes tree figures in `~/figures/introgression_figures`
