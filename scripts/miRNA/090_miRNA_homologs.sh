#!/bin/bash

makeblastdb -in db/mature.fa-no-tca.fa_dna.fa -dbtype nucl || exit 1

mkdir -p data/090_miRNA_homologs/ || exit 1
 

./091_blast_qcov_short.pl -query data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa -db db/mature.fa-no-tca.fa_dna.fa -out data/090_miRNA_homologs/tca_mature_vs_other_miRs.blastout -num_threads 100 || exit 1
