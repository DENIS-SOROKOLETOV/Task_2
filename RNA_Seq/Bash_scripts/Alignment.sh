#!/bin/bash

BASE="/mnt/f/Data"
REF="/mnt/d/BioInf/Task_2/RNA_Seq/Reference_genome"
OUTPUT="/mnt/d/BioInf/Task_2/RNA_Seq/Alignment_results"

for suffix in 71 72 73 74 75 76 83 84 85 86 87 88; do
  s="SRR104086${suffix}"
  mkdir -p "${OUTPUT}/${s}"
  echo "Выравнивание файлов эксперимента ${s}"
  hisat2 -p 6 \
  -x "${REF}/index_hisat2/fgram_v3" \
  --known-splicesite-infile "${REF}/splicesites.txt" \
  -1 "${BASE}/${s}/${s}_1_edited.fastq" \
  -2 "${BASE}/${s}/${s}_2_edited.fastq" \
  -S "${OUTPUT}/${s}/${s}_aligned.sam" \
  2> "${OUTPUT}/${s}/hisat2.log"
  echo "Перевод SAM в BAM ${s}"
  samtools view -@ 6 -bS "${OUTPUT}/${s}/${s}_aligned.sam" > "${OUTPUT}/${s}/${s}_aligned.bam"
  echo "Сортировка BAM ${s}"
  samtools sort -@ 6 "${OUTPUT}/${s}/${s}_aligned.bam" -o "${OUTPUT}/${s}/${s}_aligned_sorted.bam"
  echo "Индексирование BAM ${s}"
  samtools index "${OUTPUT}/${s}/${s}_aligned_sorted.bam"
  echo "QC по результатам выравнивания ${s}"
  samtools flagstat "${OUTPUT}/${s}/${s}_aligned_sorted.bam" > "${OUTPUT}/${s}/flagstat.txt"
done