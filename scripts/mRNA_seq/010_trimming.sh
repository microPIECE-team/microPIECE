#!/bin/bash
command -v cutadapt >/dev/null 2>&1 || { echo "I require cutadapt but it's not installed. Aborting." >&2; exit 1; }


mkdir -p data/010_trimmed_RNA
            
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/oocytes_RNAseq_rep1_1_trim.fastq.gz data/000_mRNA_raw/oocytes_RNAseq_rep1_1.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/oocytes_RNAseq_rep1_2_trim.fastq.gz data/000_mRNA_raw/oocytes_RNAseq_rep1_2.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/oocytes_RNAseq_rep2_1_trim.fastq.gz data/000_mRNA_raw/oocytes_RNAseq_rep2_1.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/oocytes_RNAseq_rep2_2_trim.fastq.gz data/000_mRNA_raw/oocytes_RNAseq_rep2_2.fastq.gz

cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/Fifth-instar_S4_L001_R1_001_trim.fastq.gz data/000_mRNA_raw/Fifth-instar_S4_L001_R1_001.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/Fifth-instar_S4_L002_R1_001_trim.fastq.gz data/000_mRNA_raw/Fifth-instar_S4_L002_R1_001.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/Fifth-instar_S4_L003_R1_001_trim.fastq.gz data/000_mRNA_raw/Fifth-instar_S4_L003_R1_001.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/Fifth-instar_S4_L004_R1_001_trim.fastq.gz data/000_mRNA_raw/Fifth-instar_S4_L004_R1_001.fastq.gz
 
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/last-larval-instar_S1_L001_R1_001_trim.fastq.gz data/000_mRNA_raw/last-larval-instar_S1_L001_R1_001.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/last-larval-instar_S1_L002_R1_001_trim.fastq.gz data/000_mRNA_raw/last-larval-instar_S1_L002_R1_001.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/last-larval-instar_S1_L003_R1_001_trim.fastq.gz data/000_mRNA_raw/last-larval-instar_S1_L003_R1_001.fastq.gz
cutadapt -a AGATCGGAAGAGC -o data/010_trimmed_RNA/last-larval-instar_S1_L004_R1_001_trim.fastq.gz data/000_mRNA_raw/last-larval-instar_S1_L004_R1_001.fastq.gz

