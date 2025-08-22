#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_tw_lineage_vcf.%A_%a.oe
#SBATCH --job-name="CbTwLineageVcf"
#SBATCH --array=0-2



source activate CT_PopGen

Cb_VCF_raw="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250626.hard-filter.isotype.vcf.gz"
lineage_info="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_by_lineage/Taiwan_lineage/Taiwan_lineage.csv"

output_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_by_lineage/Taiwan_lineage"
mkdir -p "$output_dir"
cd "$output_dir"

categories=("TD1" "Tropical" "TH")
category="${categories[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date +"%F %T")] Processing category: $category"

bcftools view \
    -S <( awk -F',' -v cat="$category" '$3==cat{print $1}' "$lineage_info" ) \
    -Oz "$Cb_VCF_raw" \
    -o "${category}.vcf.gz"


bcftools index "${category}.vcf.gz"

echo "[$(date +"%F %T")] Done $category"


