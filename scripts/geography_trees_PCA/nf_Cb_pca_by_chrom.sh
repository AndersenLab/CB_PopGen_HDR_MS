

########################################################################
############################## Chr I ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/


tmux new -s PCA_I

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/
mkdir -p I
cd I


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs/I.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/I

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA_I










########################################################################
############################## Chr II ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/


tmux new -s PCA_II

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/
mkdir -p II
cd II


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs/II.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/II

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA_II













########################################################################
############################## Chr III ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/


tmux new -s PCA_III

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/
mkdir -p III
cd III


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs/III.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/III

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA_III













########################################################################
############################## Chr IV ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/


tmux new -s PCA_IV

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/
mkdir -p IV
cd IV


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs/IV.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/IV

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA_IV











########################################################################
############################## Chr V ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/


tmux new -s PCA_V

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/
mkdir -p V
cd V


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs/V.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/V

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA_V










########################################################################
############################## Chr X ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/


tmux new -s PCA_X

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/
mkdir -p X
cd X


##### upload the eigpar_no_removal and eigpar_outliers_removed file 
#### beacuse smartpca excludes the X chromosome by default.


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs/X.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/Cb_pruned_VCF_and_PCA/sample_list.txt \
--species c_briggsae \
--vcf_folder $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/Nikita_PopGen_Brig_Project/2025_PopGen_Bri/processed_data/PCA_by_chrom/X \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA_X












