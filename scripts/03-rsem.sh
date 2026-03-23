#!/bin/bash

module add rsem

REF=ref

for R1 in fastq/*R1.fastq.gz; do
SAMPLE="$(basename $R1 | sed -e 's/_R1.*//')"
echo "$SAMPLE (rsem)"
mkdir -p align/$SAMPLE

if [ -e align/$SAMPLE/$SAMPLE.Aligned.toTranscriptome.out.bam ]; then

mkdir -p align/$SAMPLE/rsem

cgsub -t 12:00:00 -p 8 -m 32G - \
rsem-calculate-expression \
  --alignments \
  --paired-end \
  -p 8 \
  align/$SAMPLE/$SAMPLE.Aligned.toTranscriptome.out.bam \
  ref/rsem/GRCh38p14 \
  align/$SAMPLE/rsem

fi

done
