#!/bin/bash

IN_DIR="/home/rhooads/sra/aligned"
OUT_DIR="/home/rhooads/sra/counts"
INDEX_DIR="/home/rhooads/genomes"
ANNO="Homo_sapiens.GRCh38.83.gtf"

cd $IN_DIR

for sample in *sam
do
describer=$(echo ${sample} | sed 's/_Aligned.out.sam//')
echo $sample
echo $describer

echo $describer start_HTSeq-Counts

htseq-count --stranded=no $IN_DIR/${sample} $INDEX_DIR/$ANNO > $OUT_DIR/${describer}.cnt

echo $describer end_HTSeq-Counts

done
