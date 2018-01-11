#!/bin/bash
mkdir -p data/021_completed_mature_miRBase/
./021_parse_miRBaseOUT.pl data/012_miRDeep2_output_bwt1/result-bwt1.csv db/tca_mature_mirbase.fa tca-mir-3811c-1,tca-mir-3811c-2,tca-mir-3851a-1,tca-mir-3851a-2 > data/021_completed_mature_miRBase/tca_mature_mirbase_completed.fa
