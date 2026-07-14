#!/usr/bin/env bash
set -euo pipefail

for dataset in all hdr; do
    src_root="../phylo_py_pruner/pruner_${dataset}"
    dest_root="mast_${dataset}"

    mkdir -p "$dest_root"

    for rundir in "${src_root}"/occ*_smintaxa_1to1prune_quartet_*_mins75; do
        [[ -d "$rundir" ]] || continue

        matrix="${rundir}/phylopypruner_output/supermatrix.fas"
        [[ -f "$matrix" ]] || continue

        base=$(basename "$rundir")
        # base = occ80_JU3200_smintaxa_1to1prune_quartet_01_mins75

        prefix=${base%%_smintaxa_*}
        # prefix = occ80_JU3200

        quartet=$(echo "$base" | sed -E 's/.*_(quartet_[0-9]+)_mins75/\1/')
        # quartet = quartet_01

        ln -s "${PWD}/${matrix}" \
            "${dest_root}/${prefix}_${quartet}_supermatrix.fas"
    done
done
