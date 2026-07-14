#activate environment
source activate of3_env

orthofinder -M msa -S diamond_ultra_sens -A mafft -n $prefix -t 48 -f $prot_dir

