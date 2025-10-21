#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=Cb_vcf_to_zarr_geo.oe
#SBATCH --job-name="CbV2ZaGeo"


# source activate vcf_zarr
source activate CT_PopGen

#Define the inputs and outputs
## use full paths

out_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data"
vcf_raw_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/geo_vcf"
script_path="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/pi_theta_d_by_geo"



### copy the latest version scripts
cp -rn $HOME/vast-eande106/projects/Bowen/software/pi_theta_d_python_v20250530/* $script_path


### Copy and Index the VCF ###
#make the output directories
mkdir -p $out_dir/pi_theta_d_geo
mkdir -p $out_dir/pi_theta_d_geo/vcf_and_zarr

# get the file name of the VCF
for vcf_file in "$vcf_raw_dir"/*.vcf.gz; do
    vcf_name=$(basename "$vcf_file")
    ln -s "$vcf_file" "$out_dir/pi_theta_d_geo/vcf_and_zarr/$vcf_name"
done


# get the index of the VCF
for index_file in "$vcf_raw_dir"/*.vcf.gz.tbi; do
    index_name=$(basename "$index_file")
    cp -n "$index_file" "$out_dir/pi_theta_d_geo/vcf_and_zarr/$index_name"
done


# # index the VCF
# for vcf_file in "$out_dir"/pi_theta_d_geo/vcf_and_zarr/*.vcf.gz; do
#     bcftools index -f "$vcf_file"
# done




source activate vcf_zarr

vcf_files=($out_dir/pi_theta_d_geo/vcf_and_zarr/*.vcf.gz)

# for loop
for vcf_file in "${vcf_files[@]}"; do
    # extract file names
    vcf_name=$(basename "$vcf_file")

    # zarr file path
    zarr="$out_dir/pi_theta_d_geo/vcf_and_zarr/${vcf_name%.vcf.gz}.zarr"

    # export following message
    echo "How many CPUs you asked for ${SLURM_NPROCS}"
    echo "Processing file: $vcf_file"
    echo "Output Zarr file: $zarr"

    # run python script
    python $script_path/parallel_testing.py \
        --nproc "${SLURM_NPROCS}" \
        --vcf "$vcf_file" \
        --out "$zarr"

    echo "Finished processing $vcf_file"
done
