#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=Calculate_Cb_run_admixture_arrays_%A_%a.oe
#SBATCH --job-name="CbAdmixture"
#SBATCH --array=1-290


cd ../../processed_data
mkdir -p Cb_admixture
cd Cb_admixture

out_folder="$PWD/"

ln -s "../../processed_data/Cb_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/INPUTFILES/eigenstrat_input.ped" "${out_folder}LD_0.9.ped"

seed_file="${out_folder}seeds.txt"
if [ -f "$seed_file" ]; then
    echo "Loading seeds from $seed_file..."
    readarray -t seeds < "$seed_file"
else
    echo "Error: $seed_file not found."
    exit 1
fi

pops=($(seq 2 30))
total_pops=${#pops[@]}
total_seeds=${#seeds[@]}

expected_tasks=$((total_pops * total_seeds))
if [ $SLURM_ARRAY_TASK_COUNT -ne $expected_tasks ]; then
    echo "Error: Expected $expected_tasks tasks but SLURM_ARRAY_TASK_COUNT is $SLURM_ARRAY_TASK_COUNT"
    exit 1
fi

task_id=$SLURM_ARRAY_TASK_ID
task_idx=$((task_id - 1))

pop_idx=$((task_idx / total_seeds))
seed_idx=$((task_idx % total_seeds))

if [ $pop_idx -ge $total_pops ]; then
    echo "Error: pop_idx $pop_idx out of range (0-$((total_pops-1)))"
    exit 1
fi

current_pop=${pops[$pop_idx]}
current_seed=${seeds[$seed_idx]}

echo "Processing K=${current_pop} with seed=${current_seed} (task $task_id of $SLURM_ARRAY_TASK_COUNT)"


task_dir="${out_folder}K${current_pop}_seed${current_seed}"
mkdir -p "$task_dir"
cd "$task_dir"

if [ ! -f "input.ped" ]; then
    ln -s "${out_folder}LD_0.9.ped" "input.ped"
fi


source activate ADMIXTURE

echo "Running ADMIXTURE for K=${current_pop} with seed=${current_seed}"
admixture --cv=10 -s "$current_seed" "input.ped" $current_pop -j24 | tee "log_${current_pop}_${current_seed}.out"

if [ $? -eq 0 ]; then
    cp "input.${current_pop}.P" "${out_folder}LD_0.9_${current_pop}_${current_seed}.P"
    cp "input.${current_pop}.Q" "${out_folder}LD_0.9_${current_pop}_${current_seed}.Q"
    echo "Successfully completed K=${current_pop} seed=${current_seed}"
else
    echo "ADMIXTURE failed for K=${current_pop} seed=${current_seed}"
    exit 1
fi
