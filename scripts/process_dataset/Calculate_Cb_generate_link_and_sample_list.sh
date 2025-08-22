#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 00:01:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_generate_link_and_sample_list.oe
#SBATCH --job-name="CbLinList"



### generate link
ln -s $HOME/vast-eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF
ln -s $HOME/vast-eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz.tbi $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF

source activate CT_PopGen

raw_VCF="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250331.hard-filter.isotype.vcf.gz"


#####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/
mkdir -p PCA_all
cd PCA_all

bcftools query -l $raw_VCF > sample_list.txt


