#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_subset_raw_vcf.oe
#SBATCH --job-name="CbSbVCF"



source activate CT_PopGen

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/

VCF_IN="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/Archive/Archive_20250829/WI.20250626.hard-filter.isotype.vcf.gz"
VCF_OUT="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"

### remove 4 isotypes and monoallelic variants

bcftools view -s "^MY681,ECA1146,JU356,ECA1503" ${VCF_IN} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${VCF_OUT}

bcftools index -t ${VCF_OUT}



