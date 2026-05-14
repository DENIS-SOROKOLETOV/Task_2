#!/bin/bash

BASE="/mnt/f/Data"
REF="/mnt/d/BioInf/Task_2/RNA_Seq/Reference_genome"
INPUT="/mnt/d/BioInf/Task_2/RNA_Seq/Alignment_results"
OUTPUT="/mnt/d/BioInf/Task_2/RNA_Seq/Quantification_results"

for suffix in 71 72 73 74 75 76 83 84 85 86 87 88; do
  s="SRR104086${suffix}"
  mkdir -p "${OUTPUT}/${s}"
  echo "Квантификация файлов эксперимента ${s}"
  featureCounts -T 6 -p -t exon -g gene_id -s 2 \
  -a "${REF}/Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.62.gtf" \
  -o "${OUTPUT}/${s}/counts.txt" "${INPUT}/${s}/${s}_aligned_sorted.bam"
done