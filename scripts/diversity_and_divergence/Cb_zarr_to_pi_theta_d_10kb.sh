#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Cb_zarr_to_pi_theta_d.oe
#SBATCH --job-name="CbZ2PTD"


cd ../../processed_data/
mkdir -p pi_theta_d
cd pi_theta_d


out_dir="../../processed_data/pi_theta_d"
chrom_length="../../data/briggsae_genome_files/05.07.22_cb_chrom_dict"
rep_bool="../../data/Cb_repeat_region.pkl"
zarr="../../processed_data/pi_theta_d/zarr/WI.20250626.hard_filter.715_isotype.vcf.gz.zarr"
script_path="../../scripts/diversity_and_divergence/pi_theta_d_python"

source activate vcf_stats


python $script_path/pi_theta_d.py \
--zarr $zarr \
--chrom_lengths $chrom_length \
--mask_missing T \
--mask_repeats T \
--repeats_file $rep_bool \
--out $out_dir \
--filter_monomorphic T

