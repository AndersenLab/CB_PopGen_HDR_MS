#!/bin/bash
#SBATCH --job-name="TCbGTR09"
#SBATCH --output=test_GTR_LD09_model.oe
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48




cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p test_GTR_LD09
cd test_GTR_LD09

cp -n $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/test_GTR_LD09

input_file="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/test_GTR_LD09/phy_file_LD_0.9.phy"


# run iqtree 
source activate tree_v2.4
iqtree -s $input_file -mset GTR -bb 1000 -T 48



