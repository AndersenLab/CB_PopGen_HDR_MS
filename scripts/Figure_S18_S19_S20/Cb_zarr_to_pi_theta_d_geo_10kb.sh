#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=Cb_zarr_to_pi_theta_d_geo.oe
#SBATCH --job-name="CbZ2Pitd"




# make new dir
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p pi_theta_d_geo
cd pi_theta_d_geo


#Define the inputs and outputs
out_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_geo"
chrom_length="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/05.07.22_cb_chrom_dict"
rep_bool="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/Cb_repeat_region.pkl"
zarr="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d_geo/vcf_and_zarr"
script_path="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/pi_theta_d_by_geo"

# activate conda env
source activate vcf_stats


cd $zarr
mv Central\ America.zarr Central_America.zarr
mv South\ America.zarr South_America.zarr
mv North\ America.zarr North_America.zarr


zarr_dirs=$(find . -type d -name "*.zarr")

for zarr_dir_name in $zarr_dirs; do
    dir_name=$(basename "$zarr_dir_name"| sed 's/\.zarr$//')
    out_dir_geo="$out_dir/$dir_name"
    mkdir -p "$out_dir_geo"
    python "$script_path/pi_theta_d.py" --zarr "$zarr_dir_name" --chrom_lengths "$chrom_length" --mask_missing T --mask_repeats T --repeats_file $rep_bool --filter_monomorphic T --out "$out_dir_geo"
done


