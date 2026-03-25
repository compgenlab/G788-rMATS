#!/bin/bash

REF=ref

for R1 in fastq/*R1.fastq.gz; do
SAMPLE="$(basename $R1 | sed -e 's/_R1.*//')"
R2="$(echo $R1 | sed -e 's/_R1/_R2/')"

echo "$SAMPLE $R1/$R2 (STAR)"
mkdir -p align/$SAMPLE

if [ ! -e align/$SAMPLE/$SAMPLE.sortedByCoord.out.bam ]; then

cgsub -t 12:00:00 -p 8 -m 48G \
STAR \
  --runMode alignReads \
  --genomeDir $REF \
  --readFilesIn $R1 $R2 \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix align/$SAMPLE/$SAMPLE. \
  --twopassMode Basic \
  --outSAMattributes NH HI AS NM MD \
  --outSAMunmapped Within \
  --quantMode TranscriptomeSAM \
  --runThreadN 8 \; samtools index align/$SAMPLE/$SAMPLE.sortedByCoord.out.bam

fi

done
