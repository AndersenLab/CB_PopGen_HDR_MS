#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Calculate_Dxy_generate_VCF.oe
#SBATCH --job-name="CbDxyVcf"


source activate CT_PopGen

cd ../../processed_data
mkdir -p Dxy_Tropical
cd Dxy_Tropical

mkdir -p vcf


Tropical_AD_list="../../processed_data/Dxy_Tropical/Tropical_AD.txt"
Tropical_KD_list="../../processed_data/Dxy_Tropical/Tropical_KD.txt"
Tropical_TD1_list="../../processed_data/Dxy_Tropical/Tropical_TD1.txt"
Tropical_Temperate_list="../../processed_data/Dxy_Tropical/Tropical_Temperate.txt"
Tropical_TH_list="../../processed_data/Dxy_Tropical/Tropical_TH.txt"

vcf_raw="../../data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"

vcf_Tro_AD="./vcf/WI.20250626.hard_filter.715_isotype_Tro_AD.vcf.gz"
vcf_Tro_KD="./vcf/WI.20250626.hard_filter.715_isotype_Tro_KD.vcf.gz"
vcf_Tro_TD1="./vcf/WI.20250626.hard_filter.715_isotype_Tro_TD1.vcf.gz"
vcf_Tro_Temperate="./vcf/WI.20250626.hard_filter.715_isotype_Tro_Temperate.vcf.gz"
vcf_Tro_TH="./vcf/WI.20250626.hard_filter.715_isotype_Tro_TH.vcf.gz"


bcftools view -S <(tail -n +1 "${Tropical_AD_list}" | cut -f1) ${vcf_raw} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${vcf_Tro_AD}
bcftools index -t ${vcf_Tro_AD}


bcftools view -S <(tail -n +1 "${Tropical_KD_list}" | cut -f1) ${vcf_raw} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${vcf_Tro_KD}
bcftools index -t ${vcf_Tro_KD}


bcftools view -S <(tail -n +1 "${Tropical_TD1_list}" | cut -f1) ${vcf_raw} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${vcf_Tro_TD1}
bcftools index -t ${vcf_Tro_TD1}


bcftools view -S <(tail -n +1 "${Tropical_Temperate_list}" | cut -f1) ${vcf_raw} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${vcf_Tro_Temperate}
bcftools index -t ${vcf_Tro_Temperate}


bcftools view -S <(tail -n +1 "${Tropical_TH_list}" | cut -f1) ${vcf_raw} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${vcf_Tro_TH}
bcftools index -t ${vcf_Tro_TH}





