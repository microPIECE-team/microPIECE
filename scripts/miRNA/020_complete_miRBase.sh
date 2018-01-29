#!/bin/bash
mkdir -p data/021_completed_mature_miRBase/ || exit 1
./021_parse_miRDeep2_output.pl -mirdeep_out data/012_miRDeep2_output_bwt1/result-bwt1.csv -mature_fasta db/tca_mature_mirbase.fa -precursor_copies tca-mir-3811c-1,tca-mir-3811c-2,tca-mir-3851a-1,tca-mir-3851a-2 > data/021_completed_mature_miRBase/tca_mature_mirbase_completed.fa || exit 1
