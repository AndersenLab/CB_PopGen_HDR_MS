#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_admixture_post_processing.oe
#SBATCH --job-name="CbAdmPostPro"

cd ../../processed_data
mkdir -p genetic_similarity_and_admixture
cd genetic_similarity_and_admixture

out_folder="../../processed_data/genetic_similarity_and_admixture/"

echo "Starting post-processing..."
grep -h CV ../../processed_data/genetic_similarity_and_admixture/K*/log_*.out | \
cut -f3- -d" " | \
sed 's/[(|):]//g' | \
sort -k1n | \
awk 'BEGIN{OFS="\t"; print "K", "CV"}; {print $0}' > admix_replicates_CV.tsv

echo "finish post-processing"


