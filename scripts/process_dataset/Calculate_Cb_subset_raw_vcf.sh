#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Cb_subset_raw_vcf.oe
#SBATCH --job-name="CbSbVCF"

source activate CT_PopGen
cd ../../data/VCF/

### WI.20250626.hard-filter.isotype.vcf.gz is publicly available in CaeNDR
### https://caendr.org/data/data-release/c-briggsae/latest
VCF_IN="../../data/VCF/WI.20250626.hard-filter.isotype.vcf.gz"
VCF_OUT="../../data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz"

### remove 4 isotypes and monoallelic variants
bcftools view -s "^MY681,ECA1146,JU356,ECA1503" ${VCF_IN} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${VCF_OUT}

bcftools index -t ${VCF_OUT}

