#!/bin/bash
#SBATCH -A mschatz1
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=pal2nal_array
#SBATCH --output=logs/pal2nal_%A_%a.out
#SBATCH --error=logs/pal2nal_%A_%a.err

source activate pal2nal

msa=$(sed -n "${SLURM_ARRAY_TASK_ID}p" msa_files.txt)

base=$(basename "$msa" .aligned.fa)

cds="../cds_fasta/${base}.cds.fa"
out="../codon_aware_msa_wgaps/${base}.aligned.codons.gapped.fa"

if [[ ! -f "$cds" ]]; then
  echo "Missing CDS file: $cds" >&2
  exit 1
fi

pal2nal.pl "$msa" "$cds" \
  -output fasta \
  > "$out"
