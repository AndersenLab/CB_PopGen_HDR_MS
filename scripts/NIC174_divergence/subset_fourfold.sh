#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=ff_subset.oe
#SBATCH --job-name="ff_subset"

source activate bcftools

bcftools view \
  -R fourfold_closed_intervals.3col \
  -Oz \
  -o strain.MY754_JU2220.fourfold.vcf.gz \
  --threads 48 \
  strain.MY754_JU2220.vcf.gz
