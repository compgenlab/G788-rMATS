#!/bin/bash

cd $(dirname $0)

echo "file	pair	repl	sample	tissue	url" > pairs.txt
tabl export 'File accession','Paired end','Biological replicate(s)','Experiment accession','Biosample term name','File download URL' samples.txt | tail -n +2 >> pairs.txt

while read -r file pair repl sample tissue url; do
if [ $file = "file" ]; then
	continue
fi
FQ=$file.fastq.gz
FQ2="fastq/${sample}_rep${repl}_R${pair}.fastq.gz"
if [ -e "$FQ2" ]; then
	echo "$FQ2 downloaded"
	continue
fi
if [ -e "$FQ" ]; then
	mv $FQ $FQ2
	continue
fi
curl -LO $url > $FQ2.tmp && mv $FQ2.tmp $FQ2
done < pairs.txt


cd ref/
if [ ! -e GRCh38.p14.genome.fa.gz ]; then
curl -LO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.p14.genome.fa.gz
fi
cd ..


LIVER=""
SKIN=""
while read -r file pair repl sample tissue url; do
if [ "$pair" = "1" ]; then
BAM="align/${sample}_rep${repl}/${sample}_rep${repl}.Aligned.sortedByCoord.out.bam"
LIVER="$LIVER,$BAM"
fi
done < <(grep 'liver' pairs.txt)

while read -r file pair repl sample tissue url; do
if [ "$pair" = "1" ]; then
BAM="align/${sample}_rep${repl}/${sample}_rep${repl}.Aligned.sortedByCoord.out.bam"
SKIN="$SKIN,$BAM"
fi
done < <(grep 'skin' pairs.txt)

LIVER="${LIVER#,}"
SKIN="${SKIN#,}"

echo $LIVER > liver.txt
echo $SKIN > skin.txt

