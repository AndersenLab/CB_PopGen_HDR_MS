#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_make_geo_vcf_non_cosmopolitan.oe
#SBATCH --job-name="CbVcfNonCosm"



source activate CT_PopGen

Cb_VCF_raw="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250626.hard-filter.isotype.vcf.gz"
Cb_cosmopolitan_list="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Geo_info/Cb_cosmopolitan_isotype.txt"


cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p geo_vcf
cd geo_vcf


### make non_cosmopolitan vcf
bcftools view -S ^$Cb_cosmopolitan_list -Oz $Cb_VCF_raw -o Cb_non_cosmopolitan.vcf.gz
bcftools index Cb_non_cosmopolitan.vcf.gz


