#!/bin/bash

# Get the header plus your top hits by |ΔPSI|
head -1 rmats_output/SE.MATS.JCEC.txt > SE_interesting.txt

# Append top 10 most significant high-delta events
awk -F'\t' 'NR>1 && $20 < 0.05 && ($22 > 0.3 || $22 < -0.3)' \
  rmats_output/SE.MATS.JCEC.txt | \
  sort -t$'\t' -k22 -g | \
  head -10 >> SE_interesting.txt

rmats2sashimiplot \
  --b1 skin.txt \
  --b2 liver.txt \
  --event-type SE \
  -e SE_interesting.txt \
  --l1 Skin \
  --l2 Liver \
  -o sashimi_out/