#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 00:01:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_generate_link_and_sample_list.oe
#SBATCH --job-name="CbLinList"


source activate CT_PopGen

raw_VCF="../../data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"

cd ../../processed_data/
mkdir -p Cb_pruned_VCF_and_PCA
cd Cb_pruned_VCF_and_PCA

bcftools query -l $raw_VCF > sample_list.txt
