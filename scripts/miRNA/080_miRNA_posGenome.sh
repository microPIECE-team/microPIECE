#!/bin/bash

makeblastdb -in db/GCF_000002335.3_Tcas5.2_genomic.fna -dbtype nucl

mkdir -p data/080_miRNA_pos/

./081_blast_qcov.pl -p blastn -query db/tca_precursor_mirbase_completed_novel.fa -db db/GCF_000002335.3_Tcas5.2_genomic.fna -out data/080_miRNA_pos/ -num_threads 100





