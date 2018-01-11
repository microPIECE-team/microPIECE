#!/bin/bash
mkdir -p data/041_miRDeep_completed_with_novels/ 

./041_curated_mirdeep2fasta.pl data/012_miRDeep2_output_bwt1/result-bwt1.csv 10

# create miRNA mature dataset : miRBase + missing mature arms + novels
cat data/021_completed_mature_miRBase/tca_mature_mirbase_completed.fa result-bwt1-mature.fa > data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa
mv result-bwt1-mature.fa data/041_miRDeep_completed_with_novels/
# create miRNA precursir dataset: miRBase + novels
cat db/tca_precursor_mirbase.fa  result-bwt1-hairpin.fa > data/041_miRDeep_completed_with_novels/tca_precursor_mirbase_completed_novel.fa
mv result-bwt1-hairpin.fa data/041_miRDeep_completed_with_novels/
