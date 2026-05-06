source activate bcftools

bcftools view -s $strain $vcf |\
bcftools filter -i 'GT="alt"' -Oz -o $strain.vcf.gz
bedtools coverage -a $windows -b $strain.vcf.gz -counts > $strain.variant_counts.tsv
