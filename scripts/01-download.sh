#!/bin/bash

cd $(dirname $0)
mkdir -p fastq

echo "tissue	file	pair	repl	sample	url" > pairs.txt
while read -r line; do
        echo "heart     $line" >> pairs.txt
done < <(grep 'File accession\|heart\|atrium' samples.txt | tabl export 'File accession','Paired end','Biological replicate(s)','Experiment accession','File download URL' - | tail -n +2)
while read -r line; do
        echo "skin      $line" >> pairs.txt
done < <(grep 'File accession\|skin' samples.txt | tabl export 'File accession','Paired end','Biological replicate(s)','Experiment accession','File download URL' - | tail -n +2)

while read -r tissue file pair repl sample url; do
FQ2="fastq/${tissue}_${sample}_rep${repl}_R${pair}.fastq.gz"
if [ -e "$FQ2" ]; then
        echo "$FQ2 downloaded"
        continue
fi
curl --silent -L $url > $FQ2.tmp && mv $FQ2.tmp $FQ2 &
done < <(tail -n +2 pairs.txt)


HEART=""
SKIN=""
while read -r tissue file pair repl sample url; do
if [ "$pair" = "1" ]; then
BAM="align/${tissue}_${sample}_rep${repl}/${tissue}_${sample}_rep${repl}.Aligned.sortedByCoord.out.bam"
HEART="$HEART,$BAM"
fi
done < <(grep 'heart' pairs.txt)

while read -r tissue file pair repl sample url; do
if [ "$pair" = "1" ]; then
BAM="align/${tissue}_${sample}_rep${repl}/${tissue}_${sample}_rep${repl}.Aligned.sortedByCoord.out.bam"
SKIN="$SKIN,$BAM"
fi
done < <(grep 'skin' pairs.txt)

HEART="${HEART#,}"
SKIN="${SKIN#,}"

echo $HEART > heart.txt
echo $SKIN > skin.txt