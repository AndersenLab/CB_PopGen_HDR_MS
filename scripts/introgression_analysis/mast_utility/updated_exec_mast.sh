#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <mast_script>"
    exit 1
fi

mast_script="$1"

topodir="quartet_topologies"
matrixdirs=("mast_all" "mast_hdr")
#matrixdirs=("mast_hdr")

for matrixdir in "${matrixdirs[@]}"; do
    setname=$(basename "$matrixdir")

    for supermat in "$matrixdir"/occ*_quartet_*_supermatrix.clean.fas; do
        [[ -e "$supermat" ]] || continue

        base=$(basename "$supermat")

        occ=$(echo "$base" | grep -o 'occ[0-9]\+')
        tag=$(echo "$base" | grep -oE 'JU3200|JU3237|QG4232')
        quartet=$(echo "$base" | grep -o 'quartet_[0-9]\+')

        topo=$(ls "$topodir"/${quartet}_*_${tag}_topologies.tre 2>/dev/null | head -n 1)

        if [[ -z "$topo" ]]; then
            echo "Missing topology for ${quartet} ${tag}" >&2
            continue
        fi

        prefix="${setname}_${occ}_${tag}_${quartet}"

        sbatch \
          --export="supermat=$(realpath "$supermat"),topo=$(realpath "$topo"),prefix=$prefix" \
          "$mast_script"
    done
done
