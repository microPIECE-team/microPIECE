#! /bin/bash
command -v cutadapt >/dev/null 2>&1 || { echo "I require cutadapt but it's not installed. Aborting." >&2; exit 1; }

cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o 
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o 
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o 

cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o 
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o 
cutadapt -a GTGTCAGTCACTTCCAGCGG -m 20 --trim-n -o 
