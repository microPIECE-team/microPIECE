#!/bin/bash

#get maximum threads
#NPROC="$(nproc)"

mkdir -p data/020_mapped_RNA

gunzip db/GCF_000002335.3_Tcas5.2_genomic.fna.gz
gunzip db/GCF_000002335.3_Tcas5.2_genomic.gff.gz
# build genome index
/opt/hisat2-2.1.0/hisat2-build --seed 1337 db/GCF_000002335.3_Tcas5.2_genomic.fna db/TCA52_HISAT2

# align raw reads to genome
/opt/hisat2-2.1.0/hisat2 --seed 1337 -p 100 --dta-cufflinks -x db/TCA52_HISAT2 -1 data/000_mRNA_raw/oocytes_RNAseq_rep1_1.fastq.gz,data/000_mRNA_raw/oocytes_RNAseq_rep2_1.fastq.gz -2 data/000_mRNA_raw/oocytes_RNAseq_rep1_2.fastq.gz,data/000_mRNA_raw/oocytes_RNAseq_rep2_2.fastq.gz -S data/020_mapped_RNA/oocyte.sam

/opt/hisat2-2.1.0/hisat2 --seed 1337 -p 100 --dta-cufflinks -x db/TCA52_HISAT2 -U data/000_mRNA_raw/Fifth-instar_S4_L001_R1_001.fastq.gz,data/000_mRNA_raw/Fifth-instar_S4_L002_R1_001.fastq.gz,data/000_mRNA_raw/Fifth-instar_S4_L003_R1_001.fastq.gz,data/000_mRNA_raw/Fifth-instar_S4_L004_R1_001.fastq.gz -S data/020_mapped_RNA/L5.sam

/opt/hisat2-2.1.0/hisat2 --seed 1337 -p 100 --dta-cufflinks -x db/TCA52_HISAT2 -U data/000_mRNA_raw/last-larval-instar_S1_L001_R1_001.fastq.gz,data/000_mRNA_raw/last-larval-instar_S1_L002_R1_001.fastq.gz,data/000_mRNA_raw/last-larval-instar_S1_L003_R1_001.fastq.gz,data/000_mRNA_raw/last-larval-instar_S1_L004_R1_001.fastq.gz -S data/020_mapped/LLI.sam

samtools view -bS -@ 100 data/020_mapped_RNA/oocyte.sam | samtools sort -@ 100 - -o data/020_mapped_RNA/oocyte.bam
samtools view -bS -@ 100 data/020_mapped_RNA/L5.sam | samtools sort -@ 100 - -o data/020_mapped_RNA/L5.bam
samtools view -bS -@ 100 data/020_mapped/LLI.sam | samtools sort -@ 100 - -o data/020_mapped/LLI.bam

