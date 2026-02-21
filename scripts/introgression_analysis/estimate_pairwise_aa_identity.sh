#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name="run_pimm"

if [[ -z "$input" || -z "$OUTDIR" ]]; then
    echo "ERROR: Must provide input and OUTDIR via --export"
    echo "Usage: sbatch --export=input=file.fa,OUTDIR=/path/to/out estimate_pairwise_aa_identity.sh"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$OUTDIR"

base=$(basename "$input" .fa)
outfile="$OUTDIR/${base}.pimm"

echo "Processing $input"
echo "Output → $outfile"

awk -v outfile="$outfile" '
BEGIN {
    OFS=","
}
/^>/ {
    # store previous seq if any
    if (seq != "") { seqs[id]=seq; seq="" }
    id=substr($0,2)
    next
}
{
    seq=seq $0
}
END {
    # store last sequence
    if (seq != "") seqs[id]=seq

    # collect ids
    n=0
    for (i in seqs) { ids[n]=i; n++ }

    # header row
    printf "," > outfile
    for (i=0; i<n; i++) printf ids[i] (i<n-1?",":"\n") >> outfile

    # pairwise identity (symmetric counting)
    for (i=0; i<n; i++) {
        s1 = seqs[ids[i]]
        printf ids[i]"," >> outfile
        for (j=0; j<n; j++) {
            s2 = seqs[ids[j]]

            l1 = length(s1)
            l2 = length(s2)
            L  = (l1 > l2 ? l1 : l2)

            matches = 0
            valid   = 0

            for (k = 1; k <= L; k++) {
                a = (k <= l1 ? substr(s1,k,1) : "-")
                b = (k <= l2 ? substr(s2,k,1) : "-")

                # skip positions where either is a gap
                if (a == "-" || b == "-") continue

                valid++
                if (a == b) matches++
            }

            iden = (valid > 0) ? matches / valid : 0
            printf iden (j < n-1 ? "," : "\n") >> outfile
        }
    }
}
' "$input"
