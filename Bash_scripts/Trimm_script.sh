#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/f/Data"

for suffix in 71 72 73 74 75 76; do
  s="SRR104086${suffix}"
  mkdir -p "${BASE}/${s}"
  echo "Обработка файлов эксперимента ${s}"
  fastp \
    --detect_adapter_for_pe \
    -i "${BASE}/${s}/${s}_1.fastq" -I "${BASE}/${s}/${s}_2.fastq" \
    -o "${BASE}/${s}/${s}_1_wo_adapters.fastq" -O "${BASE}/${s}/${s}_2_wo_adapters.fastq"
  fastp \
    -i "${BASE}/${s}/${s}_1_wo_adapters.fastq" -I "${BASE}/${s}/${s}_2_wo_adapters.fastq" \
    -o "${BASE}/${s}/${s}_1_edited.fastq" -O "${BASE}/${s}/${s}_2_edited.fastq" \
    --trim_front1 12 \
    --trim_front2 12
done

for suffix in 84 85 86 87 88; do
  s="SRR104086${suffix}"
  mkdir -p "${BASE}/${s}"
  echo "Обработка файлов эксперимента ${s}"
  fastp \
    --detect_adapter_for_pe \
    -i "${BASE}/${s}/${s}_cleaned.1.fastq" -I "${BASE}/${s}/${s}_cleaned.2.fastq" \
    -o "${BASE}/${s}/${s}_1_wo_adapters.fastq" -O "${BASE}/${s}/${s}_2_wo_adapters.fastq"
  fastp \
    -i "${BASE}/${s}/${s}_1_wo_adapters.fastq" -I "${BASE}/${s}/${s}_2_wo_adapters.fastq" \
    -o "${BASE}/${s}/${s}_1_edited.fastq" -O "${BASE}/${s}/${s}_2_edited.fastq" \
    --trim_front1 12 \
    --trim_front2 12
done
