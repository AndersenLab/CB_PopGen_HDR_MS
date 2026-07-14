#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=dp_subset.oe
#SBATCH --job-name="dp_subset"

source activate bcftools

bcftools view \
    --threads 24 \
    -s MY754,JU2220 \
    -Oz \
    -o strain.MY754_JU2220.vcf.gz \
    strain.hard.vcf.gz


