#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=Cb_vcf_to_zarr.oe
#SBATCH --job-name="CbV2Zarr"


source activate CT_PopGen

cd ../../processed_data/
mkdir -p pi_theta_d
cd pi_theta_d

out_dir="../../processed_data/pi_theta_d"
raw_VCF="../../data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"

mkdir -p $out_dir/vcf
mkdir -p $out_dir/zarr

vcf_name=$(basename $raw_VCF)
ln -s $raw_VCF $out_dir/vcf/$vcf_name

bcftools index $out_dir/vcf/$vcf_name
vcf_input="../../processed_data/pi_theta_d/vcf/WI.20250626.hard_filter.715_isotype.vcf.gz"

source activate vcf_zarr

zarr=$out_dir/zarr/$vcf_name.zarr

echo How many CPUs you asked for ${SLURM_NPROCS}
python ../../scripts/diversity_and_divergence/pi_theta_d_python/parallel_testing.py \
--nproc ${SLURM_NPROCS} \
--vcf $vcf_input \
--out $zarr
