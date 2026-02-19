#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_make_geo_vcf.%A_%a.oe
#SBATCH --job-name="CbGeoVCF"
#SBATCH --array=0-11


source activate CT_PopGen

Cb_VCF_raw="../../data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"
geo_info="../../processed_data/Geo_info/Cb_indep_isotype_info_geo.csv"

output_dir="../../processed_data/geo_vcf"
mkdir -p "$output_dir"
cd "$output_dir"

categories=("Africa" "Asia" "Taiwan" "Hawaii" "Pacific" "Caribbean" "Cosmopolitan" "North America" "South America" "Central America" "Australia" "Europe")
category="${categories[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date +"%F %T")] Processing category: $category"
bcftools view \
    -S <( awk -F',' -v cat="$category" '$4==cat{print $1}' "$geo_info" ) \
    -Oz "$Cb_VCF_raw" \
    -o "${category}.vcf.gz"

bcftools index "${category}.vcf.gz"
echo "[$(date +"%F %T")] Done $category"


