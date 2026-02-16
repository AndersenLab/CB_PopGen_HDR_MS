#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Cb_iterations_5_10_15_20_25_30_35_40.oe
#SBATCH --job-name="CbIter5_40"



par_path="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/PCA_more_iterations"



###### 5 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_5_iterations
cd PCA_5_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_5_iterations.par" > smartpca_5.log

sed -n -e '/Tracy/,$p' smartpca_5.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv




###### 10 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_10_iterations
cd PCA_10_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_10_iterations.par" > smartpca_10.log

sed -n -e '/Tracy/,$p' smartpca_10.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv




###### 15 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_15_iterations
cd PCA_15_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_15_iterations.par" > smartpca_15.log

sed -n -e '/Tracy/,$p' smartpca_15.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv




###### 20 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_20_iterations
cd PCA_20_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_20_iterations.par" > smartpca_20.log

sed -n -e '/Tracy/,$p' smartpca_20.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv



###### 25 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_25_iterations
cd PCA_25_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_25_iterations.par" > smartpca_25.log

sed -n -e '/Tracy/,$p' smartpca_25.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv





###### 30 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_30_iterations
cd PCA_30_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_30_iterations.par" > smartpca_30.log

sed -n -e '/Tracy/,$p' smartpca_30.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv











###### 35 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_35_iterations
cd PCA_35_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_35_iterations.par" > smartpca_35.log

sed -n -e '/Tracy/,$p' smartpca_35.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv










###### 40 iterations #####
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_more_iterations
cd PCA_more_iterations
mkdir -p PCA_40_iterations
cd PCA_40_iterations

source activate eigensoft
smartpca -p "$par_path/Cb_par_40_iterations.par" > smartpca_40.log

sed -n -e '/Tracy/,$p' smartpca_40.log | \
  sed -e '/kurt/,$d' | \
  awk '$0 !~ /##/ && $0 !~ /^[[:space:]]*#/ {print}' | \
  sed -e 's/[[:space:]]\+/ /g' | \
  sed 's/^ //g' | \
  awk 'BEGIN {print "N\teigenvalue\tdifference\ttwstat\tp-value\teffect.n"} {print}' | \
  awk -F" " '{$1=$1; print}' OFS="\t" > TracyWidom_statistics_outlier_removal.tsv

