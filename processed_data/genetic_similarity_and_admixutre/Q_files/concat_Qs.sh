awk '
BEGIN {
    FS = "[[:space:]]+"   # input fields can be any whitespace
    OFS = "\t"            # output fields separated by tabs
}

# First file: sample IDs
NR == FNR {
    samples[++n] = $0
    next
}

# All Q files
{
    # Get clean filename and group
    fname = FILENAME
    sub(/^.*\//, "", fname)   # strip path
    sub(/\.Q$/, "", fname)    # strip extension
    split(fname, parts, "_")

    # K = third underscore field
    k = parts[3] + 0

    # Build the ancestry proportion fields, normalized to tabs
    q = ""
    for (j = 1; j <= k; j++) {
        q = q (j > 1 ? OFS : "") $j
    }

    # Write to group-specific file
    out = "concat_Qfiles_K" parts[3] ".tsv"
    print samples[FNR], q, fname >> out
}
' ../../HDRs/CB_all_strain_list.txt *.Q
