#!/usr/bin/env bash

outfile="mast_summary.tsv"

printf "set\tocc\ttag\tquartet\tsubmodel\tstatus\tBIC\tTW1\tTW2\tTW3\n" > "$outfile"

for logfile in mast_*_occ*_*_quartet_*.submodel*.log; do
    [[ -e "$logfile" ]] || continue

    base="${logfile%.log}"

    set=$(echo "$base" | grep -oE '^mast_(all|hdr)')
    occ=$(echo "$base" | grep -oE 'occ[0-9]+')
    tag=$(echo "$base" | sed -E 's/^mast_(all|hdr)_occ[0-9]+_([^_]+)_quartet_[0-9]+\.submodel[0-9]+$/\2/')
    quartet=$(echo "$base" | grep -oE 'quartet_[0-9]+')
    submodel=$(echo "$base" | grep -oE 'submodel[0-9]+')

    iqfile="${base}.iqtree"

    if grep -q "ERROR" "$logfile"; then
        status="ERROR"
    elif [[ -f "$iqfile" ]]; then
        status="OK"
    else
        status="NO_IQTREE"
    fi

    if [[ -f "$iqfile" ]]; then
        bic=$(grep "Bayesian information criterion (BIC) score:" "$iqfile" | awk '{print $NF}')

        weights=$(grep "Tree weights:" "$iqfile" | sed 's/.*Tree weights: //; s/,//g')
        read -r tw1 tw2 tw3 <<< "$weights"

        [[ -z "$bic" ]] && bic="NA"
        [[ -z "$tw1" ]] && tw1="NA"
        [[ -z "$tw2" ]] && tw2="NA"
        [[ -z "$tw3" ]] && tw3="NA"
    else
        bic="NA"
        tw1="NA"
        tw2="NA"
        tw3="NA"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$set" "$occ" "$tag" "$quartet" "$submodel" "$status" \
        "$bic" "$tw1" "$tw2" "$tw3" \
        >> "$outfile"
done
