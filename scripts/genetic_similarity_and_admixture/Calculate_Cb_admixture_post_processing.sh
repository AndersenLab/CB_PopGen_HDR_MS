#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_admixture_post_processing.oe
#SBATCH --job-name="CbAdmPostPro"







cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p Cb_admixture
cd Cb_admixture

out_folder="/home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_admixture/"

cp -n /home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_admixture/K*/*.out \
/home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_admixture/


echo "Starting post-processing..."
grep -h CV log_*.out | \
cut -f3- -d" " | \
sed 's/[(|):]//g' | \
sort -k1n | \
awk 'BEGIN{OFS="\t"; print "K", "CV"}; {print $0}' > admix_replicates_CV.tsv

echo "finish post-processing"


