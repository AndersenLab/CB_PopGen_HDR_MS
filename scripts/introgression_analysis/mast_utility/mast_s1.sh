#!/bin/bash
#SBATCH -A mschatz1
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=masts1.oe
#SBATCH --job-name="masts1"

source activate iqtree3
iqtree3 -s $supermat -m "TMIX{LG+FO+G,LG+FO+G,LG+FO+G}+T" -te $topo --prefix ${prefix}.submodel1 -T 48
