#!/bin/bash
time ./005_trimming.sh &&
time ./010_mining.sh &&
time ./020_complete_miRBase.sh &&
time ./040_mirdeep2fasta.sh &&
time ./060_quantification.sh &&
time ./070_isomiR.sh &&
time ./080_miRNA_posGenome.sh &&
time ./090_miRNA_homologs.sh 
