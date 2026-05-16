#!/usr/bin/env bash

DATA="/mnt/d/BioInf/WGBS_data"
SAVE="/mnt/d/BioInf/Task_2/WGBS/QC_results"
data_dirs=("$DATA"/*/)
# suffix="_*_edited"

for dir_name in ${data_dirs[@]}; do
	cd $dir_name
	fastq_files=("*.fq")
	mkdir -p "${SAVE}/$(basename ${dir_name})"
	fastqc ${fastq_files[0]} ${fastq_files[1]} -o "${SAVE}/$(basename ${dir_name})"
	cd ..
done