#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name="nucmer"

#Given a target genome and a query genome, generate TSV with nucmer alignments between the two genomes

#$prefix = output prefix
#$ref = path_to reference genome
#$query = path_to wild genome

#activate environment
source activate nucmer

# align with nucmer (will spit out a .delta file)
nucmer --maxgap=500 --mincluster=100 --prefix=$prefix --coords $ref $query

# get coordinate file
show-coords -r -l -T $prefix.delta | awk '$5 > 1000' > ${prefix}_transformed.tsv 
