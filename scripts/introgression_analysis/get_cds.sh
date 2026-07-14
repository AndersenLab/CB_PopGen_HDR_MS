#includes gffread
source activate gfftools

#loop through GFF, produce spliced CDS fasta
for gff in ../${dir}/*.gff
do
    sample=$(basename $gff .braker.gff3)
    genome="../genomes/${sample}.genome.fa"
    out="${sample}.spliced_cds.fa"

    if [[ -f $genome ]]; then
        echo "Processing $sample"
        gffread -x $out -g $genome $gff
    else
        echo "Genome not found for $sample: $genome" >&2
    fi
done
