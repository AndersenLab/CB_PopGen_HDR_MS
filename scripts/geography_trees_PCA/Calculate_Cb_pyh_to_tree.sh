#!/bin/bash
#SBATCH --job-name="CbPyh2Tree"
#SBATCH --output=Calculate_Cb_pyh_to_tree.%A_%a.oe
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --array=0-3




cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/LD_pruned_trees


######SBATCH --array=0-3

# define input files 
files=(
    "/home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy"
    "/home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/LD_pruned_trees/LD_0.8/phy_file_LD_0.8.phy"
    "/home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/LD_pruned_trees/LD_0.7/phy_file_LD_0.7.phy"
    "/home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/LD_pruned_trees/LD_0.6/phy_file_LD_0.6.phy"
)

# set array ID 
input_file=${files[$SLURM_ARRAY_TASK_ID]}


# run iqtree 
source activate tree_v2.4
iqtree -s $input_file -mset GTR -bb 1000 -T 48



