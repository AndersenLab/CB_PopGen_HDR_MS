#!/bin/bash

#SBATCH -J IntProScan                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 48                            # Number of cores

output="../../processed_data/gene_enrichment"
mkdir -p $output
TMP=${SLURM_TMPDIR:-/scratch4/eande106/Lance/ipr.$SLURM_JOB_ID}
mkdir -p $TMP

/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/database/interproscan-5.75-106.0/interproscan.sh \
	--formats TSV \
	--input ../../processed_data/gene_enrichment/c_briggsae.QX1410_20250929.csq.longest.nomito.prot.fa \
	--goterms \
	--cpu 48 \
	--iprlookup \
	--disable-precalc \
	--output-file-base $output/QX_IPR_allApps_20251006 \
	--tempdir $TMP
