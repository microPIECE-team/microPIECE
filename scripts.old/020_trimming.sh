#! /bin/bash
command -v cutadapt >/dev/null 2>&1 || { echo "I require cutadapt but it's not installed. Aborting." >&2; exit 1; }
mkdir ../020/
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o ../020/SRR5163635_trim.fastq.gz ../data/AAE/SRR5163635.fastq.gz
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o ../020/SRR5163636_trim.fastq.gz ../data/AAE/SRR5163636.fastq.gz
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o ../020/SRR5163637_trim.fastq.gz ../data/AAE/SRR5163637.fastq.gz

cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o ../020/SRR5163632_trim.fastq.gz ../data/AAE/SRR5163632.fastq.gz
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o ../020/SRR5163633_trim.fastq.gz ../data/AAE/SRR5163633.fastq.gz
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o ../020/SRR5163634_trim.fastq.gz ../data/AAE/SRR5163634.fastq.gz
