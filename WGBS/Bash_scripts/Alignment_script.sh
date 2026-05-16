#!/bin/bash

REF="/mnt/d/BioInf/Task_2/RNA_Seq/Reference_genome/Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.dna.toplevel.fa"
OUTPUT="/mnt/d/BioInf/WGBS/Alignment_results"
THREADS=6

for suffix in 75 76 77 78; do
  s="SRR103828${suffix}"
  mkdir -p "${OUTPUT}/${s}"
  R1=""
  R2=""
  for DATA_ROOT in "/mnt/d/BioInf/WGBS_data" "/mnt/e/BioInf/WGBS_data" "/mnt/e/WGBS_data"; do
    for ext in fq fq.gz fastq fastq.gz; do
      cand1="${DATA_ROOT}/${s}/${s}_1_val_1.${ext}"
      cand2="${DATA_ROOT}/${s}/${s}_2_val_2.${ext}"
      if [[ -f "${cand1}" && -f "${cand2}" ]]; then
        R1="${cand1}"
        R2="${cand2}"
        break 2
      fi
    done
  done
  if [[ -z "${R1}" ]]; then
    echo "Не найдены обрезанные пары ${s}_1_val_1.* / ${s}_2_val_2.* (fq|fq.gz|fastq|fastq.gz) в /mnt/d/BioInf/WGBS_data, /mnt/e/BioInf/WGBS_data, /mnt/e/WGBS_data" >&2
    exit 1
  fi
  echo "Выравнивание и сортировка BAM для эксперимента ${s}"
  bwameth.py --threads "${THREADS}" \
  --reference "${REF}" \
  --read-group "${s}" \
  "${R1}" \
  "${R2}" \
  -Y \
  2> "${OUTPUT}/${s}/bwameth.log" \
  | samtools view -@ "${THREADS}" -b \
  | samtools sort -@ "${THREADS}" -o "${OUTPUT}/${s}/${s}_aligned_sorted.bam" -
  echo "Индексирование BAM ${s}"
  samtools index "${OUTPUT}/${s}/${s}_aligned_sorted.bam"
  echo "QC по результатам выравнивания ${s}"
  samtools flagstat "${OUTPUT}/${s}/${s}_aligned_sorted.bam" > "${OUTPUT}/${s}/flagstat.txt"
done
