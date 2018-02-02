#!/bin/bash
REQ_PROGS=(cutadapt bwa samtools bedtools)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

#get maximum threads
NPROC="$(nproc)"


mkdir -p data/001_trim_smRNA/ || exit 1
#trimming the raw sequencing reads
./006_cutadapt_folder.pl data/000_raw_smRNA/ data/001_trim_smRNA/ || exit 1

#filtering the trimmed reads by ncRNAs except miRNAs
mkdir -p data/002_filter_smRNA/ || exit 1

# create database from ncRNAs
gunzip -f -k db/TCA_all_ncRNA_but_miR.fa_dna.fa.gz || exit 1
bwa index db/TCA_all_ncRNA_but_miR.fa_dna.fa || exit 1

# alignment of smallRNA to ncRNA
for i in data/001_trim_smRNA/* 
do 
bwa aln -n 1 -o 0 -e 0 -k 1 -t $NPROC -f "$i".ncrna.sai db/TCA_all_ncRNA_but_miR.fa_dna.fa "$i" || exit 1;
bwa samse -f "$i".sam db/TCA_all_ncRNA_but_miR.fa_dna.fa "$i".ncrna.sai "$i" || exit 1;
samtools view -b -f 4 "$i".sam | bedtools bamtofastq -i - -fq "$i".ncrna_ual.fq || exit 1;
rm "$i".ncrna.sai "$i".sam || exit 1;
mv "$i".ncrna_ual.fq data/002_filter_smRNA/ || exit 1;
done



