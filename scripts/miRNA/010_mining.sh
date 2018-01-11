#!/bin/bash
REQ_PROGS=(miRDeep2.pl bwa)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

gunzip -f -k db/GCF_000002335.3_Tcas5.2_genomic.fna.gz
gunzip -f -k db/tca_mature_mirbase.fa.gz
gunzip -f -k db/mature.fa-no-tca.fa.gz
gunzip -f -k db/tca_precursor_mirbase.fa.gz



mkdir -p data/011_concat_smRNA/
cat data/002_filter_smRNA/* > data/011_concat_smRNA/TCA_smallRNA_concat.fastq


#mkdir -p output_bwa
mkdir -p data/012_miRDeep2_output_bwt1

#./011_miRDeep2_bwa.pl data/003_concat_smRNA/ output_bwa/ db/GCF_000002335.3_Tcas5.2_genomic.fna db/tca_mature_mirbase.fa db/mature.fa-no-tca.fa db/tca_precursor_mirbase.fa 100
#mv result_*.csv result-bwa.csv

./012_miRDeep2_bwt1.pl data/011_concat_smRNA/ data/012_miRDeep2_output_bwt1/ db/GCF_000002335.3_Tcas5.2_genomic.fna db/tca_mature_mirbase.fa db/mature.fa-no-tca.fa db/tca_precursor_mirbase.fa 110
mv result_*.csv data/012_miRDeep2_output_bwt1/result-bwt1.csv

#clean up all files 
mv mapper_logs data/012_miRDeep2_output_bwt1
mv mirdeep_runs data/012_miRDeep2_output_bwt1
mv bowtie.log data/012_miRDeep2_output_bwt1
mv expression* data/012_miRDeep2_output_bwt1
mv miRNAs_expressed_all_samples_* data/012_miRDeep2_output_bwt1
mv pdfs_* data/012_miRDeep2_output_bwt1
mv dir_prepare* data/012_miRDeep2_output_bwt1
mv error* data/012_miRDeep2_output_bwt1
mv mirna_results* data/012_miRDeep2_output_bwt1
mv result_* data/012_miRDeep2_output_bwt1
