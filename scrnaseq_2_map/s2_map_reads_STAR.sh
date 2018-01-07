#!/bin/bash

#cd into ~/sra dir containing fastq files

OUT_DIR="/home/rhooads/sra/aligned/"

INDEX_DIR="/home/rhooads/STAR_genomes"

THREADS=8

#RNA-Seq otfen uses paired reads such as the following examples in our $IN_DIR:
#ls $IN_DIR:
#FlowCell1_L003_R1.fastq
#FlowCell1_L003_R2.fastq
#FlowCell2_L007_R1.fastq
#FlowCell2_L007_R2.fastq


#select file
for sample in *_1.fastq
do
describer=$(echo ${sample} | sed 's/_..fastq//')
echo $sample
echo $describer

#select analysis
echo $describer start_STAR

#if only running non-paired RNA-Seq delete: "${describer}_R2.fastq" below
STAR --runThreadN $THREADS --genomeDir $INDEX_DIR --outFileNamePrefix $OUT_DIR --readFilesIn ${describer}_R1.fastq ${describer}_R2.fastq

echo $describer end_STAR

done