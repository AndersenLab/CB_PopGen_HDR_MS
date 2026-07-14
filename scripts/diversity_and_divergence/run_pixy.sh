#!/bin/bash
#SBATCH -A mschatz1_bigmem
#SBATCH -p bigmem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=pixy_%x_%j.oe
#SBATCH --job-name="pixy"

source activate pixy_env

OUTDIR="pixy_${PREFIX}"

mkdir -p "$OUTDIR"

pixy \
  --stats pi dxy fst watterson_theta tajima_d \
  --vcf "$VCF" \
  --populations populations.txt \
  --window_size 10000 \
  --n_cores 48 \
  --output_folder "$OUTDIR" \
  --output_prefix "$PREFIX"
