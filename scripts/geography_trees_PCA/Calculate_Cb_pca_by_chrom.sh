#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --output=Calculate_Cb_vcf_by_chrom.%a.oe
#SBATCH --job-name="CbVcfChrom"
#SBATCH --array=0-5


source activate CT_PopGen

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/
mkdir -p PCA_by_chrom
cd PCA_by_chrom


input_vcf="$HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"


CHROMS=(I II III IV V X)





### if no id, use first Chrom
IDX=${SLURM_ARRAY_TASK_ID:-0}
CHR=${CHROMS[$IDX]}

OUTDIR="./by_chrom_vcfs"
mkdir -p "${OUTDIR}"
OUT_VCF="${OUTDIR}/${CHR}.vcf.gz"
LOG="${OUTDIR}/${CHR}.bcftools.log"


echo "Running: bcftools view --regions ${CHR} ${input_vcf} -Oz -o ${OUT_VCF}"
bcftools view --regions "${CHR}" "${input_vcf}" -Oz -o "${OUT_VCF}" 2> "${LOG}"
tabix -p vcf "${OUT_VCF}"



echo "Finished chromosome ${CHR} at $(date)"
echo "Output: ${OUT_VCF} (index ${OUT_VCF}.tbi)"

