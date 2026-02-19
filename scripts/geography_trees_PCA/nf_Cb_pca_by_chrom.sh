# This script need to be run in a cluster with a SLURM scheduler
# Here is a link to the pipeline in use for further detail:
# https://github.com/AndersenLab/post-gatk-nf

cd ../../processed_data/PCA_by_chrom/
mkdir -p ${chr}
cd ${chr}

nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf ../../processed_data/PCA_by_chrom/by_chrom_vcfs/X.vcf.gz \
--pops ../../processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder ../../processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output ../../processed_data/PCA_by_chrom/${chr} \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal

#### The following parameters in the eigpar_no_removal and eigpar_outliers_removed file are required
#### --noxdata: NO 
#### Beacuse smartpca excludes the X chromosome by default.

