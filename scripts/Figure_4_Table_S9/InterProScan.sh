#!/bin/bash

#SBATCH -J IntProScan                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 48                            # Number of cores

wkdir="../../processed_data/gene_enrichment"

./interproscan.sh \
	--formats TSV \
	--input $wkdir/c_briggsae.QX1410_20250929.csq.longest.nomito.prot.fa \
	--goterms \
	--cpu 48 \
	--iprlookup \
	--disable-precalc \
	--output-file-base $wkdir/QX_IPR_allApps_20251006 \
	--tempdir $TMP
