#!/bin/bash

# run gsnap to map the reads of CLIP sequencing to AAE genome

# check for required programs
REQ_PROGS=(gsnap gmap_build samtools bedtools)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

# creation of temporary folder for this step
mkdir ../030/

# build genome index for gsnap
gmap_build -g -D ../030/gsnap_db -k 15 -d Aedes_aegypti ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz

# gsnap mapping including conversion to bam to skip temporary sam-file creation
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163635_trim.fastq.gz | samtools view -bS - | samtools sort - ../030/SRR5163635_trim_gsnap.bam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163636_trim.fastq.gz | samtools view -bS - | samtools sort - ../030/SRR5163636_trim_gsnap.bam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163637_trim.fastq.gz | samtools view -bS - | samtools sort - ../030/SRR5163637_trim_gsnap.bam

gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163632_trim.fastq.gz | samtools view -bS - | samtools sort - ../030/SRR5163632_trim_gsnap.bam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163633_trim.fastq.gz | samtools view -bS - | samtools sort - ../030/SRR5163633_trim_gsnap.bam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163634_trim.fastq.gz | samtools view -bS - | samtools sort - ../030/SRR5163634_trim_gsnap.bam

# bam2bed
bedtools bamtobed -i ../030/SRR5163635_trim_gsnap.bam > ../030/SRR5163635_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163636_trim_gsnap.bam > ../030/SRR5163636_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163637_trim_gsnap.bam > ../030/SRR5163637_trim_gsnap.bed

bedtools bamtobed -i ../030/SRR5163632_trim_gsnap.bam > ../030/SRR5163632_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163633_trim_gsnap.bam > ../030/SRR5163633_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163634_trim_gsnap.bam > ../030/SRR5163634_trim_gsnap.bed


