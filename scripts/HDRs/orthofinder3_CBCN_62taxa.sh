#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --job-name="of3_ft"


#given the path to a directory with protein FASTA files, run orthofinder and retrieve orthogroups file (Orthogroups.tsv)

#activate environment
source activate of3_env

orthofinder -M msa -S diamond_ultra_sens -A mafft -n CBCN_mafft_ft -t 48 -f ../../processed_data/gene_diversity/prot_CBCN_62_taxa/

#mv ../../processed_data/gene_diversity/prot_CBCN_62_taxa/OrthoFinder/Results_CBCN_mafft_ft/Orthogroups/Orthogroups.tsv ../../processed_data/gene_diversity/CBCN_orthogroups.tsv
