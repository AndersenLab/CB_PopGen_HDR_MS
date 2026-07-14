#!/usr/bin/env bash
set -euo pipefail

src="${PWD}/paired_files"

for tag in JU3200 JU3237 QG4232; do
    list="hdr_ogs_${tag}.txt"
    dest="paired_files_hdrs_${tag}"

    mkdir -p "$dest"

    while read -r og; do
    	[[ -z "$og" ]] && continue

    	fa="${src}/${og}.fa"
	tree="${src}/${og}.treefile"

    	# Skip OGs missing either file
    	if [[ ! -L "$fa" || ! -f "$tree" ]]; then
        	echo "Skipping ${og}: missing .fa or .treefile"
        	continue
    	fi

    	cp -P "$fa" "$dest/"
	ln -s "$tree" "${dest}/${og}.treefile"

	done < "$list"
done
