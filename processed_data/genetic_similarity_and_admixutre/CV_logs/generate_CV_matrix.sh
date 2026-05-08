gawk '
BEGIN { OFS = "\t" }

# When a new file starts, parse K and seed from filename: log_<K>_<seed>.out
FNR == 1 {
    file = FILENAME
    sub(/^.*\//, "", file)              # strip path
    if (match(file, /^log_([0-9]+)_([0-9]+)\.out$/, m)) {
        K    = m[1] + 0
        seed = m[2] + 0
    } else {
        next
    }
}

# Find the CV error line and capture the numeric value
/^CV[[:space:]]*error/ {
    if (match($0, /CV[[:space:]]*error[[:space:]]*\(K=[0-9]+\):[[:space:]]*([0-9.eE+-]+)/, mm)) {
        val = mm[1] + 0
        seenK[K] = 1
        seenSeed[seed] = 1
        cv[K, seed] = val
    }
    next
}

END {
    # Collect and sort K and seed lists numerically
    nk = 0; for (k in seenK)    Ks[++nk] = k
    ns = 0; for (s in seenSeed) Seeds[++ns] = s
    asort(Ks, Ks, "@val_num_asc")
    asort(Seeds, Seeds, "@val_num_asc")

    # Header
    printf "K"
    for (i = 1; i <= ns; i++) printf OFS "S" Seeds[i]
    printf "\n"

    # Rows
    for (i = 1; i <= nk; i++) {
        k = Ks[i]
        printf "%s", k
        for (j = 1; j <= ns; j++) {
            s = Seeds[j]
            if ((k, s) in cv) printf OFS "%s", cv[k, s]
            else              printf OFS ""
        }
        printf "\n"
    }
}
' log_*.out > ../cv_matrix.tsv
