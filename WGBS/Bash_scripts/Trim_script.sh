#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/e/WGBS_data"
OUTPUT="/mnt/d/BioInf/WGBS_data"

for suffix in 75 76 77 78; do
  s="SRR103828${suffix}"
  mkdir -p "${OUTPUT}/${s}"
  echo "Обработка файлов эксперимента ${s}"
  trim_galore --paired \
    --illumina \
    --stringency 5 \
    --quality 20 \
    --length 30 \
    --output_dir "${OUTPUT}/${s}" \
    "${BASE}/${s}/${s}_1.fastq" "${BASE}/${s}/${s}_2.fastq"
done
