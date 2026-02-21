#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --job-name="of3_ft"


#given the path to a directory with protein FASTA files, run orthofinder and retrieve orthogroups file (Orthogroups.tsv), resolved gene trees (Resolved_Gene_Trees.txt), and consensus tree (consensus_species_tree.txt).

#activate environment
source activate of3_env

orthofinder -M msa -S diamond_ultra_sens -A famsa -n CBCN_famsa_trop_ft -t 48 -f ../../processed_data/tree_topology/prot_TropicalCB_CN_36taxa

#mv ../../processed_data/tree_topology/prot_TropicalCB_CN_36taxa/OrthoFinder/Results_CBCN_famsa_trop_ft/Orthogroups/Orthogroups.tsv ../../processed_data/tree_topology/orthogroups_CBCN_Tropical.tsv
#mv ../../processed_data/tree_topology/prot_TropicalCB_CN_36taxa/OrthoFinder/Results_CBCN_famsa_trop_ft/Resolved_Gene_Trees/Resolved_Gene_Trees.txt ../../processed_data/tree_topology/Resolved_Gene_Trees.txt
#mv ../../processed_data/tree_topology/prot_TropicalCB_CN_36taxa/OrthoFinder/Results_CBCN_famsa_trop_ft/Species_Tree/SpeciesTree_rooted_node_labels.txt ../../processed_data/tree_topology/consensus_species_tree.txt

