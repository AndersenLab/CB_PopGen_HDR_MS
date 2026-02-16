#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Cb_zarr_to_pi_theta_d.oe
#SBATCH --job-name="CbZ2PTD"



cp -n $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/make_Cb_repeats_bed_file/Cb_repeat_region.pkl $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data
cp -n $HOME/vast-eande106/projects/Ryan/cb_pop_gen/scripts/inputs/05.07.22_cb_chrom_dict $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data


cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/
mkdir -p pi_theta_d
cd pi_theta_d


out_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d"
chrom_length="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/05.07.22_cb_chrom_dict"
rep_bool="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/Cb_repeat_region.pkl"
zarr="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d/zarr/WI.20250626.hard_filter.715_isotype.vcf.gz.zarr"
script_path="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/pi_theta_d"

source activate vcf_stats


python $script_path/pi_theta_d.py \
--zarr $zarr \
--chrom_lengths $chrom_length \
--mask_missing T \
--mask_repeats T \
--repeats_file $rep_bool \
--out $out_dir \
--filter_monomorphic T

