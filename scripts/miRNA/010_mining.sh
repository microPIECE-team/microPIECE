#!/bin/bash
REQ_PROGS=(miRDeep2.pl bwa)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

gunzip -f -k db/GCF_000002335.3_Tcas5.2_genomic.fna.gz || exit 1
gunzip -f -k db/tca_mature_mirbase.fa.gz || exit 1
gunzip -f -k db/mature.fa-no-tca.fa.gz || exit 1
gunzip -f -k db/tca_precursor_mirbase.fa.gz || exit 1



mkdir -p data/011_concat_smRNA/
cat data/002_filter_smRNA/* > data/011_concat_smRNA/TCA_smallRNA_concat.fastq


mkdir -p data/012_miRDeep2_output_bwt1


./012_miRDeep2_bwt1.pl -dir data/011_concat_smRNA/ -out data/012_miRDeep2_output_bwt1/ -ref_genome db/GCF_000002335.3_Tcas5.2_genomic.fna -species_mature_miRs db/tca_mature_mirbase.fa -other_mature_miRs db/mature.fa-no-tca.fa -species_precursor_mirs db/tca_precursor_mirbase.fa -threads 110
mv result_*.csv data/012_miRDeep2_output_bwt1/result-bwt1.csv

#clean up all files 
mv mapper_logs mirdeep_runs bowtie.log expression* miRNAs_expressed_all_samples_* pdfs_* dir_prepare* error* mirna_results* result_* -t data/012_miRDeep2_output_bwt1
