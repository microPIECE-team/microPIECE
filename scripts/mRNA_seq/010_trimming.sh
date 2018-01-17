#!/bin/bash
command -v cutadapt >/dev/null 2>&1 || { echo "I require cutadapt but it's not installed. Aborting." >&2; exit 1; }


mkdir -p data/010_trimmed_RNA

for i in data/000_mRNA_raw/*.fastq.gz
do
   echo "Working on file '${i}'"
   
   FILENAME=$(basename "${i}" .fastq.gz)
   
   cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/"${FILENAME}"_trim.fastq.gz data/000_mRNA_raw/"${FILENAME}".fastq.gz
done

