#!/bin/bash
#SBATCH --job-name="CbVCF2Pyh"
#SBATCH --output=Calculate_Cb_vcf_to_pyh_%A_%a.out
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --array=0-3

cd ../../processed_data/
mkdir -p LD_pruned_trees
cd LD_pruned_trees

# define input files 
files=(
    "../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/INPUTFILES/eiganstrat_input.vcf.gz"
    "../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.8/INPUTFILES/eiganstrat_input.vcf.gz"
    "../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.7/INPUTFILES/eiganstrat_input.vcf.gz"
    "../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.6/INPUTFILES/eiganstrat_input.vcf.gz"
)

input_file=${files[$SLURM_ARRAY_TASK_ID]}
ld_value=$(basename $(dirname $(dirname $input_file)))

mkdir -p $ld_value
cd $ld_value

python ../../scripts/vcf2phylip.py -i $input_file
mv eiganstrat_input.min4.phy phy_file_${ld_value}.phy
