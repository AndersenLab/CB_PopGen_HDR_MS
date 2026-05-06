source activate bcftools

### WI.20260324.hard-filter.isotype.vcf.gz is publicly available in CaeNDR
### https://caendr.org/data/data-release/c-briggsae/latest
VCF_IN="../../data/VCF/WI.20260324.hard-filter.isotype.vcf.gz"
VCF_OUT="../../data/VCF/WI.20260324.hard-filter.713_isotype.vcf.gz"

### remove 4 isotypes and monoallelic variants
bcftools view -s "^MY681,ECA1146,JU356,ECA1503" ${VCF_IN} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' -Oz -o ${VCF_OUT}

bcftools index -t ${VCF_OUT}

