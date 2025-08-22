#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=Cb_vcf_to_zarr.oe
#SBATCH --job-name="CbV2Zarr"


source activate CT_PopGen

# Creat a new dir
cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/
mkdir -p pi_theta_d
cd pi_theta_d

### copy the latest version scripts
cp -rn $HOME/vast-eande106/projects/Bowen/software/pi_theta_d_python_v20250530/* $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/pi_theta_d/



#Define the inputs and outputs
## use full paths
out_dir="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d"
raw_VCF="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250626.hard-filter.isotype.vcf.gz"


### Copy and Index the VCF ###
#make the output directories
mkdir -p $out_dir/vcf
mkdir -p $out_dir/zarr


# get the file name of the VCF
vcf_name=$(basename $raw_VCF)
ln -s $raw_VCF $out_dir/vcf/$vcf_name


# index the VCF
bcftools index $out_dir/vcf/$vcf_name

vcf_input="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/pi_theta_d/vcf/WI.20250626.hard-filter.isotype.vcf.gz"



source activate vcf_zarr

#define the zarr output path relative to the conversion directory
zarr=$out_dir/zarr/$vcf_name.zarr



echo How many CPUs you asked for ${SLURM_NPROCS}
python /home/bwang97/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/scripts/pi_theta_d/parallel_testing.py \
--nproc ${SLURM_NPROCS} \
--vcf $vcf_input \
--out $zarr


