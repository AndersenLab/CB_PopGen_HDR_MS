#!/bin/bash
#SBATCH -A mschatz1
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=masts6.oe
#SBATCH --job-name="masts6"

source activate iqtree3
iqtree3 -s $supermat -m "LG+FO+G+T" -te $topo --prefix ${prefix}.submodel6 -T 48
