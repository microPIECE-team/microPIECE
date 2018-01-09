#!/bin/bash

./041_curated_mirdeep2fasta.pl result-bwt1.csv 10

# create miRNA mature dataset : miRBase + missing mature arms + novels
cat tca_mature_mirbase_completed.fa result-bwt1-mature.fa > db/tca_mature_mirbase_completed_novel.fa
rm -rf tca_mature_mirbase_completed.fa
# create miRNA precursir dataset: miRBase + novels
cat db/tca_precursor_mirbase.fa  result-bwt1-hairpin.fa > db/tca_precursor_mirbase_completed_novel.fa

