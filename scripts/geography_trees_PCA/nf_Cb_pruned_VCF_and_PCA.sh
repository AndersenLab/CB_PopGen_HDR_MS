########################################################################
############################## PCA all ###################################
########################################################################

# This script need to be run in a cluster with a SLURM scheduler
# Here is a link to the pipeline in use for further detail:
# https://github.com/AndersenLab/post-gatk-nf

source activate /data/eande106/software/conda_envs/nf24_env

cd ../../processed_data
mkdir -p Cb_pruned_VCF_and_PCA
cd Cb_pruned_VCF_and_PCA

nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf ../../data/VCF/WI.20250626.hard_filter.715_isotype.vcf.gz \
--pops ../../processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder ../../data/VCF \
--eigen_ld 0.9,0.8,0.7,0.6 \
--postgatk false \
--singletons false \
--delly false \
--output ../../processed_data/Cb_pruned_VCF_and_PCA \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal

#### The following parameters in the eigpar_no_removal and eigpar_outliers_removed file are required
#### --noxdata: NO 
#### Beacuse smartpca excludes the X chromosome by default.
