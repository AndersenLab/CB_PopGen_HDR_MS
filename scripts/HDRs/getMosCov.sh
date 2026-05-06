source activate bcftools

mosdepth $strain $bam_dir/${strain}.bam -b $windows -t 5 -T 1,2,5 -n
