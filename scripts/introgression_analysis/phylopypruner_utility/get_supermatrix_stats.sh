#!/usr/bin/env bash
set -euo pipefail

out="supermatrix_stats_summary.tsv"

printf "Gene set\tCB Background\tQuartet ID\tNumber of alignments\tPercent missing data\tSupermatrix length\n" > "$out"

for geneset in ALL HDRs; do
    if [[ "$geneset" == "ALL" ]]; then
        root="pruner_all"
    else
        root="pruner_hdr"
    fi

    for stats in "$root"/occ*_smintaxa_1to1prune_quartet_*_mins75/phylopypruner_output/supermatrix_stats.csv; do
        [[ -f "$stats" ]] || continue

        rundir=$(basename "$(dirname "$(dirname "$stats")")")

        cb=$(echo "$rundir" | grep -oE 'JU3200|JU3237|QG4232')
        quartet=$(echo "$rundir" | grep -oE 'quartet_[0-9]+')

        awk -F',' -v geneset="$geneset" -v cb="$cb" -v quartet="$quartet" '
            $1 == "output" {
                print geneset "\t" cb "\t" quartet "\t" $2 "\t" $10 "\t" $11
            }
        ' "$stats" >> "$out"
    done
done
