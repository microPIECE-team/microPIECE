#!/bin/bash
REQ_PROGS=(miRDeep2.pl bwa RNAfold)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done


zcat ../../data/TCA/GCF_000002335.3_Tcas5.2_genomic.fna.gz > db/GCF_000002335.3_Tcas5.2_genomic.fna

mkdir output_bwa
mkdir output_bwt1

./miRDeep2_bwa.pl input/ output_bwa/ db/GCF_000002335.3_Tcas5.2_genomic.fna db/tca_mature_mirbase.fa db/mature.fa-no-tca.fa db/tca_precursor_mirbase.fa 20


./miRDeep2_bwt1.pl input/ output_bwt1/ db/GCF_000002335.3_Tcas5.2_genomic.fna db/tca_mature_mirbase.fa db/mature.fa-no-tca.fa db/tca_precursor_mirbase.fa 20
