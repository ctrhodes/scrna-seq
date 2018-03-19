#!/bin/bash

IN_DIR="/home/rhooads/sra"
OUT_DIR="/home/rhooads/sra/aligned"
INDEX_DIR="/home/rhooads/genomes/h38_ercc.idx"

RD1_SUFFIX="_1"
RD2_SUFFIX="_2"

cd $IN_DIR


#fastq or fastq.gz? Update for header, sed and input suffixes
for sample in *_1.fastq.gz
do
describer=$(echo ${sample} | sed 's/_..fastq.gz//')
echo $sample
echo $describer

echo $describer start_kallisto-quant

kallisto quant -i $INDEX_DIR -o $OUT_DIR/${describer} -b 50 ${describer}$RD1_SUFFIX.fastq.gz ${describer}$RD2_SUFFIX.fastq.gz

echo $describer end_kallisto-quant

done
