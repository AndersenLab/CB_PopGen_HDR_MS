#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=Cb_zarr_to_pi_theta_d_tw_lineage.oe
#SBATCH --job-name="CbZ2PitdTwLineage"




# make new dir
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_by_lineage
mkdir -p Taiwan_lineage
cd Taiwan_lineage


#Define the inputs and outputs
out_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_by_lineage/Taiwan_lineage"
chrom_length="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/05.07.22_cb_chrom_dict"
rep_bool="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/Cb_repeat_region.pkl"
zarr="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_by_lineage/Taiwan_lineage/vcf_and_zarr"
script_path="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/isotype_by_lineage"

# activate conda env
source activate vcf_stats

### raw script ###
# python scripts/scripts/pi_theta_d.py --zarr $zarr --chrom_lengths $chrom_length --mask_missing T --out $out_dir

## make it a for-loop
# find all dir end with .zarr
cd $zarr


zarr_dirs=$(find . -type d -name "*.zarr")

# loop every dir
for zarr_dir_name in $zarr_dirs; do
    # extract dir name (.zarr)
    dir_name=$(basename "$zarr_dir_name"| sed 's/\.zarr$//')
    # define export dir
    out_dir_geo="$out_dir/$dir_name"
    # Create output directory if it doesn't exist
    mkdir -p "$out_dir_geo"
    # run the script
    python "$script_path/pi_theta_d.py" --zarr "$zarr_dir_name" --chrom_lengths "$chrom_length" --mask_missing T --mask_repeats T --repeats_file $rep_bool --filter_monomorphic T --out "$out_dir_geo" --window_size 1000
done


