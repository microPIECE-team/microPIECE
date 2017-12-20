#!/bin/bash
mkdir -p data/060_miRNA_expression/
./061_sam2de.pl 062_s2d_cfg db/tca_mature_mirbase_completed_novel.fa > data/060_miRNA_expression/TCA_miRNA_expression.csv
