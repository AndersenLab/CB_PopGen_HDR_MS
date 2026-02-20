#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --output=Calculate_Dxy_Tropical_%A_%a.oe
#SBATCH --job-name="CbDxyTropical"
#SBATCH --array=0-4 


source activate pixy

cd ../../processed_data/diversity_and_divergence/
mkdir -p Dxy_Tropical
cd Dxy_Tropical

VCF_DIR=../../processed_data/diversity_and_divergence/Dxy_Tropical/vcf
POP_DIR=../../processed_data/diversity_and_divergence/Dxy_Tropical

VCF_FILES=(
    "WI.20250626.hard_filter.715_isotype_Tro_AD.vcf.gz"
    "WI.20250626.hard_filter.715_isotype_Tro_KD.vcf.gz"
    "WI.20250626.hard_filter.715_isotype_Tro_TD1.vcf.gz"
    "WI.20250626.hard_filter.715_isotype_Tro_Temperate.vcf.gz"
    "WI.20250626.hard_filter.715_isotype_Tro_TH.vcf.gz"
)

POP_FILES=(
    "Tropical_AD.txt"
    "Tropical_KD.txt"
    "Tropical_TD1.txt"
    "Tropical_Temperate.txt"
    "Tropical_TH.txt"
)

VCFI=${VCF_DIR}/${VCF_FILES[$SLURM_ARRAY_TASK_ID]}
POP_FILE=${POP_DIR}/${POP_FILES[$SLURM_ARRAY_TASK_ID]}

pixy --stats dxy \
  --vcf $VCFI \
  --populations $POP_FILE \
  --output_folder $POP_DIR \
  --window_size 10000 \
  --bypass_invariant_check yes \
  --output_prefix $(basename ${POP_FILE%.txt})_ \
  --n_cores 15 

