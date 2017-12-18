#!/bin/bash

mkdir -p data/020_mapped_RNA

gunzip db/GCF_000002335.3_Tcas5.2_genomic.fna.gz
gunzip db/GCF_000002335.3_Tcas5.2_genomic.gff.gz
# build genome index
bowtie2-build db/GCF_000002335.3_Tcas5.2_genomic.fna db/TCA52_BWT2


tophat2 -p 100 -o data/020_mapped_RNA/oocyte/ -G db/GCF_000002335.3_Tcas5.2_genomic.gff db/TCA52_BWT2 data/010_trimmed_RNA/oocytes_RNAseq_rep1_1_trim.fastq.gz,data/010_trimmed_RNA/oocytes_RNAseq_rep2_1_trim.fastq.gz data/010_trimmed_RNA/oocytes_RNAseq_rep1_2_trim.fastq.gz,data/010_trimmed_RNA/oocytes_RNAseq_rep2_2_trim.fastq.gz

tophat2 -p 100 -o data/020_mapped_RNA/L5/ -G db/GCF_000002335.3_Tcas5.2_genomic.gff db/TCA52_BWT2 data/010_trimmed_RNA/Fifth-instar_S4_L001_R1_001_trim.fastq.gz,data/010_trimmed_RNA/Fifth-instar_S4_L002_R1_001_trim.fastq.gz,data/010_trimmed_RNA/Fifth-instar_S4_L003_R1_001_trim.fastq.gz,data/010_trimmed_RNA/Fifth-instar_S4_L004_R1_001_trim.fastq.gz

tophat2 -p 100 -o data/020_mapped_RNA/LLI/ -G db/GCF_000002335.3_Tcas5.2_genomic.gff db/TCA52_BWT2 data/010_trimmed_RNA/last-larval-instar_S1_L001_R1_001_trim.fastq.gz,data/010_trimmed_RNA/last-larval-instar_S1_L002_R1_001_trim.fastq.gz,data/010_trimmed_RNA/last-larval-instar_S1_L003_R1_001_trim.fastq.gz,data/010_trimmed_RNA/last-larval-instar_S1_L004_R1_001_trim.fastq.gz
