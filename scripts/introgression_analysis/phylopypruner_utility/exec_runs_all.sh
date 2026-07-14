#!/usr/bin/env bash
set -euo pipefail

occ="$1"

tar_dir="../paired_files/"

for exclude_file in \
  ../JU3200_quartet_exclude.tsv \
  ../JU3237_quartet_exclude.tsv \
  ../QG4232_quartet_exclude.tsv
do
    tag=$(basename "$exclude_file" _quartet_exclude.tsv)

    tail -n +2 "$exclude_file" | while IFS=$'\t' read -r quartet_id exclude exclude_otus; do
        sbatch \
          --export="quartet_id=$quartet_id,exclude_otus=$exclude_otus,occ=$occ,tar_dir=$tar_dir,tag=$tag" \
          ../../../../scripts/introgression_analysis/pruner_nopdist_smintaxa_1to1pruner_exclQ.sh
    done
done
