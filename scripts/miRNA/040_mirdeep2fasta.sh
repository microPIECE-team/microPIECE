#!/bin/bash

./041_curated_mirdeep2fasta.pl result-bwa-bwt1.csv high-conf

# create miRNA mature dataset : miRBase + missing mature arms + novels
cat tca_mature_mirbase_completed.fa result-bwa-bwt1-mature.fa > db/tca_mature_mirbase_completed_novel.fa

# create miRNA precursir dataset: miRBase + novels
cat db/tca_precursor_mirbase.fa  result-bwa-bwt1-hairpin.fa > db/tca_precursor_mirbase_completed_novel.fa
