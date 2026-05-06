#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 08:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --output=subset_vcf.oe
#SBATCH --job-name="subset_vcf"

source activate bcftools

bcftools view -S $samples -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${invcf%.713_isotype.vcf.gz}.tropical_subset.vcf.gz $invcf --threads 12
