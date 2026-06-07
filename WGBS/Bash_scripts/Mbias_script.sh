#!/bin/bash

REF="/mnt/d/BioInf/Task_2/RNA_Seq/Reference_genome/Fusarium_graminearum_ph_1_gca_000240135.ASM24013v3.dna.toplevel.fa"
INPUT="/mnt/d/BioInf/WGBS/Alignment_results"
OUTPUT="/mnt/d/BioInf/Task_2/WGBS/Extract_results/Mbias_results"
THREADS=6

if ! command -v MethylDackel >/dev/null 2>&1; then
  echo "Команда MethylDackel не найдена в PATH" >&2
  exit 1
fi

if command -v multiqc >/dev/null 2>&1; then
  MULTIQC_CMD=(multiqc)
elif command -v python3 >/dev/null 2>&1; then
  MULTIQC_CMD=(python3 -m multiqc)
elif command -v py >/dev/null 2>&1; then
  MULTIQC_CMD=(py -m multiqc)
else
  echo "Команда multiqc не найдена в PATH" >&2
  exit 1
fi

mkdir -p "${OUTPUT}"

for suffix in 75 76 77 78; do
  s="SRR103828${suffix}"
  BAM="${INPUT}/${s}/${s}_aligned_sorted.bam"
  if [[ ! -f "${BAM}" ]]; then
    echo "Не найден BAM: ${BAM}" >&2
    exit 1
  fi
  mkdir -p "${OUTPUT}/${s}"
  echo "Построение M-bias таблиц для эксперимента ${s}"
  if ! MethylDackel mbias "${REF}" "${BAM}" "${OUTPUT}/${s}/${s}" \
  --txt \
  --CHG \
  --CHH \
  -@ "${THREADS}" \
  > "${OUTPUT}/${s}/${s}_mbias.txt" 2> "${OUTPUT}/${s}/mbias.log"; then
    echo "Ошибка MethylDackel mbias для ${s}" >&2
    exit 1
  fi
done

echo "Построение сводного отчета MultiQC"
mkdir -p "${OUTPUT}/multiqc"
if ! "${MULTIQC_CMD[@]}" "${OUTPUT}" \
-o "${OUTPUT}/multiqc" \
-n "mbias_multiqc_report.html" \
> "${OUTPUT}/multiqc/multiqc.log" 2>&1; then
  echo "Ошибка построения отчета MultiQC" >&2
  exit 1
fi
