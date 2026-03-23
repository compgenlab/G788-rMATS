#!/bin/bash
READ_OVERHANG=100 # read length - 1

cd $(dirname $0)

module add star rsem

mkdir -p fastq ref/rsem align

cd ref/

if [ ! -e GRCh38.p14.genome.fa.gz ]; then
    echo "Downloading GRCh38p14"
    curl -LO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.p14.genome.fa.gz
fi
if [ ! -e gencode.v49.annotation.gtf.gz ]; then
    echo "Downloading Gencode V49 GTF"
    curl -LO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v49.annotation.gtf.gz
fi


if [ ! -e GRCh38.p14.genome.fa ]; then
    echo "Decompressing FASTA"
    zcat GRCh38.p14.genome.fa.gz > GRCh38.p14.genome.fa
fi

if [ ! -e gencode.v49.annotation.gtf ]; then
    echo "Decompressing GTF"
    zcat gencode.v49.annotation.gtf.gz > gencode.v49.annotation.gtf
fi

if [ ! -e SAindex ]; then
    echo "Creating STAR index"
    STAR \
      --runMode genomeGenerate \
      --genomeDir . \
      --genomeFastaFiles GRCh38.p14.genome.fa \
      --sjdbGTFfile gencode.v49.annotation.gtf \
      --sjdbOverhang $READ_OVERHANG \
      --runThreadN 8
fi

if [ ! -e rsem/GRCh38p14.n2g.idx.fa ]; then
    rsem-prepare-reference --gtf gencode.v49.annotation.gtf GRCh38.p14.genome.fa rsem/GRCh38p14
fi
