#!/bin/bash

REF="/mnt/d/BioInf/Task_2/RNA_Seq/Reference_genome/Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.dna.toplevel.fa"
INPUT="/mnt/d/BioInf/WGBS/Alignment_results"
OUTPUT="/mnt/d/BioInf/Task_2/WGBS/Extract_results"
THREADS=6

for suffix in 75 76 77 78; do
  s="SRR103828${suffix}"
  BAM="${INPUT}/${s}/${s}_aligned_sorted.bam"
  if [[ ! -f "${BAM}" ]]; then
    echo "Не найден BAM: ${BAM}" >&2
    exit 1
  fi
  mkdir -p "${OUTPUT}/${s}"
  echo "Извлечение метилирования для эксперимента ${s}"
  MethylDackel extract "${REF}" "${BAM}" \
  -o "${OUTPUT}/${s}/${s}" \
  --CHG \
  --CHH \
  -@ "${THREADS}" \
  > "${OUTPUT}/${s}/methyldackel.log" 2>&1
done
