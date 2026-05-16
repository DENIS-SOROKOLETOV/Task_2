#!/bin/bash
BASE="/mnt/e/WGBS_data"
SCRIPT="/mnt/d/BioInf/Task_2/WGBS/Bash_scripts"
for suffix in 78; do
  s="SRR103828${suffix}"
  cd "${BASE}/${s}"
  echo "Получение fastq файлов эксперимента ${s}"
  scratch="${BASE}/${s}/fasterq_scratch"
  mkdir -p "${scratch}"
  fasterq-dump "${s}.sra" \
  --threads 6 \
  --split-files \
  --progress \
  -O "${BASE}/${s}" \
  -t "${scratch}"
  rm -rf "${scratch}"
  cd "${SCRIPT}"
done
