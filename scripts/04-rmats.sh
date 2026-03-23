#!/bin/bash
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --time 4:00:00

mkdir -p rmats_output rmats_tmp
module add conda
conda activate rmats

rmats.py \
  --b1 skin.txt \
  --b2 liver.txt \
  --gtf ref/gencode.v49.annotation.gtf \
  -t paired \
  --readLength 100 \
  --nthread 8 \
  --od rmats_output/ \
  --tmp rmats_tmp/ \
  --libType fr-firststrand

