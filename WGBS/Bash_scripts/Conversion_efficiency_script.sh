#!/bin/bash

INPUT="/mnt/d/BioInf/Task_2/WGBS/Extract_results"
RESULT="${INPUT}/conversion_efficiency_no_spikein.txt"
status=0

echo -e "sample\tchh_weighted_methylation_percent\tconversion_efficiency_chh_percent\tcpg_weighted_methylation_percent\tconversion_efficiency_cpg_percent" > "${RESULT}"

for suffix in 75 76 77 78; do
  s="SRR103828${suffix}"
  CHH="${INPUT}/${s}/${s}_CHH.bedGraph"
  CPG="${INPUT}/${s}/${s}_CpG.bedGraph"

  if [[ ! -f "${CHH}" ]]; then
    echo "Не найден CHH файл: ${CHH}" >&2
    status=1
    continue
  fi

  if [[ ! -f "${CPG}" ]]; then
    echo "Не найден CpG файл: ${CPG}" >&2
    status=1
    continue
  fi

  row_chh=$(awk 'BEGIN{m=0;u=0} !/^track/ && NF>=6 {m+=$5;u+=$6} END{if(m+u>0){wm=(m/(m+u))*100; ce=100-wm; printf "%.6f\t%.6f", wm, ce}else{printf "NA\tNA"}}' "${CHH}")
  row_cpg=$(awk 'BEGIN{m=0;u=0} !/^track/ && NF>=6 {m+=$5;u+=$6} END{if(m+u>0){wm=(m/(m+u))*100; ce=100-wm; printf "%.6f\t%.6f", wm, ce}else{printf "NA\tNA"}}' "${CPG}")
  echo -e "${s}\t${row_chh}\t${row_cpg}" >> "${RESULT}"
done

echo "Результаты сохранены: ${RESULT}"
exit "${status}"
