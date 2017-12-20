#!/bin/bash
REQ_PROGS=(cutadapt bwa samtools bedtools)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done




mkdir -p data/001_trim_smRNA/
#trimming the raw sequencing reads
./006_cutadapt_folder.pl data/000_raw_smRNA/ data/001_trim_smRNA/

#filtering the trimmed reads by ncRNAs except miRNAs
mkdir -p data/002_filter_smRNA/

# create database from ncRNAs
gunzip -k db/TCA_all_ncRNA_but_miR.fa_dna.fa.gz
bwa index db/TCA_all_ncRNA_but_miR.fa_dna.fa

# alignment of smallRNA to ncRNA
for i in data/001_trim_smRNA/* 
do 
bwa aln -n 1 -o 0 -e 0 -k 1 -t 100 -f $i.ncrna.sai db/TCA_all_ncRNA_but_miR.fa_dna.fa $i
bwa samse -f $i.sam db/TCA_all_ncRNA_but_miR.fa_dna.fa $i.ncrna.sai $i
samtools view -b -f 4 $i.sam > $i.ncrna_ual.bam
bedtools bamtofastq -i $i.ncrna_ual.bam -fq $i.ncrna_ual.fq
rm -rf $i.ncrna.sai
rm -rf $i.sam
rm -rf $i.ncrna_ual.bam
mv $i.ncrna_ual.fq data/002_filter_smRNA/
done



