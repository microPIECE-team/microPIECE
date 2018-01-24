#!/bin/bash
mkdir -p data/040_reformat/

./041_reformat.pl data/030_cuffdiff/isoforms.fpkm_tracking db/GCF_000002335.3_Tcas5.2_genomic.gff  > data/040_reformat/isoforms.fpkm_tracking_reformat.csv
