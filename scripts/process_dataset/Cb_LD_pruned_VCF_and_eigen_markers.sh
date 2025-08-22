

########################################################################
############################## PCA all ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########



tmux new -s PCA

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data
mkdir -p PCA_all
PCA_all


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF/WI.20250626.hard-filter.isotype.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_all/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/data/VCF \
--eigen_ld 0.9,0.8,0.7,0.6 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_all \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal


#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA









