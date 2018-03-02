#! /bin/bash
# check for tools
command -v fastq-dump >/dev/null 2>&1 || { echo "I require NCBI SRA Toolkit but it's not installed. Get it on https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software ! Aborting." >&2; exit 1; }


# obtain data from public servers - may take 2-3 h

# Aedes aegypti
mkdir ../data/AAE/
# genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/015/GCF_000004015.4_AaegL3/GCF_000004015.4_AaegL3_genomic.fna.gz -O ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz
# GFF
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/015/GCF_000004015.4_AaegL3/GCF_000004015.4_AaegL3_genomic.gff.gz -O ../data/AAE/GCF_000004015.4_AaegL3_genomic.gff.gz
# Transcriptome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/015/GCF_000004015.4_AaegL3/GCF_000004015.4_AaegL3_rna.fna.gz -O ../data/AAE/GCF_000004015.4_AaegL3_rna.fna.gz
# Proteome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/015/GCF_000004015.4_AaegL3/GCF_000004015.4_AaegL3_protein.faa.gz -O ../data/AAE/GCF_000004015.4_AaegL3_protein.faa.gz
# CLIP
# 24h
fastq-dump --gzip SRR5163635 -O ../data/AAE/
fastq-dump --gzip SRR5163636 -O ../data/AAE/
fastq-dump --gzip SRR5163637 -O ../data/AAE/
# 72h
fastq-dump --gzip SRR5163632 -O ../data/AAE/
fastq-dump --gzip SRR5163633 -O ../data/AAE/
fastq-dump --gzip SRR5163634 -O ../data/AAE/

# Tribolium castaneum
mkdir ../data/TCA/
# genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.fna.gz -O ../data/TCA/GCF_000002335.3_Tcas5.2_genomic.fna.gz
# GFF
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.gff.gz -O ../data/TCA/GCF_000002335.3_Tcas5.2_genomic.gff.gz
# Transcriptome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_rna.fna.gz -O ../data/TCA/GCF_000002335.3_Tcas5.2_rna.fna.gz
# Proteome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_protein.faa.gz -O ../data/TCA/GCF_000002335.3_Tcas5.2_protein.faa.gz 
