#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=pruner_sm_1to1_noex.oe
#SBATCH --job-name="pruner_sm_1to1_noex"

source activate py_pruner_v130

gene_occ=$(awk "BEGIN {print $occ/100}")
outdir="occ${occ}_${tag}_smintaxa_1to1prune_${quartet_id}_mins75"

mkdir -p "$outdir"

phylopypruner \
  --dir $tar_dir \
  --min-len 20 \
  --prune 1to1 \
  --min-support 0.75 \
  --min-taxa 4 \
  --min-otu-occupancy 0.1 \
  --min-gene-occupancy $gene_occ \
  --outgroup PX506_longest_prot_otu_PX506 \
  --exclude $exclude_otus \
  --threads 24 \
  --output $outdir \
  --subclades ../subclades.txt \
  --no-plot \
  --overwrite
