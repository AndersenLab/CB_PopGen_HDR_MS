#!/bin/bash
#SBATCH -A eande106_bigmem
#SBATCH -p bigmem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --output="QX1410_AF16_N2.oe"
#SBATCH --job-name="briggsaeOrtho"

directory_with_proteomes=$1

orthofinder -f $directory_with_proteomes -t 48
