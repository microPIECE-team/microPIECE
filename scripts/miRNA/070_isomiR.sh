#!/bin/bash

#command -v RNAfold >/dev/null 2>&1 || { echo "I require RNAfold but it's not installed. Aborting." >&2; exit 1; }
#command -v java -jar miraligner.jar >/dev/null 2>&1 || { echo "I require miraligner.jar but it's not installed. Aborting." >&2; exit 1; }
#command -v perl -e 'use RNA::HairpinFigure qw/draw/' >/dev/null 2>&1 || { echo "I require the PERL module RNA::HairpinFigure but it's not installed. Aborting." >&2; exit 1; }

#mkdir -p data/071_filterN_smRNA/
#mkdir -p data/072_isomiR_db/

#for i in data/001_trim_smRNA/* ; 
#do 
#	FILEBASENAME=$(basename $i) 
#	./071_filter_fastq_N.pl ${i} > data/071_filterN_smRNA/${FILEBASENAME}_filterN.fastq; 
#done;
#
#
#wget -nc ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz -P db/
#gunzip -f -k db/miRNA.str.gz
#./072_create_mirbase_struct.pl -hairpin data/041_miRDeep_completed_with_novels/tca_precursor_mirbase_completed_novel.fa -mature data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa -struct db/miRNA.str
##rm tmp_hairpin.fa tmp_struct.rna
#
## copy the hairpin and structure file into database folder (necessary for miraligner)
#mv custom.str data/072_isomiR_db/miRNA.str
#cp data/041_miRDeep_completed_with_novels/tca_precursor_mirbase_completed_novel.fa data/072_isomiR_db/hairpin.fa
#

#mkdir -p data/073_isomiR_output/

# run miraligner script
#./073_seqbuster_pipe.pl -fq_path data/071_filterN_smRNA/ -db_path data/072_isomiR_db/ -out_path data/073_isomiR_output/ -species tca

# reformat output for DB 

mkdir -p data/074_isomiR_reformat/

./074_reformat_isomiRs.pl data/073_isomiR_output/egg_S8_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/egg_S8_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/egg_S8_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/egg_S8_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna egg > data/074_isomiR_reformat/egg.csv


./074_reformat_isomiRs.pl data/073_isomiR_output/First-instar-L12_S6_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L12_S6_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L12_S6_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L12_S6_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L15_S5_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L15_S5_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L15_S5_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/First-instar-L15_S5_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna L1 > data/074_isomiR_reformat/L1.csv

./074_reformat_isomiRs.pl data/073_isomiR_output/Fifth-instar_S2_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/Fifth-instar_S2_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/Fifth-instar_S2_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/Fifth-instar_S2_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna L5 > data/074_isomiR_reformat/L5.csv

./074_reformat_isomiRs.pl data/073_isomiR_output/last-larval-instar_S1_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/last-larval-instar_S1_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/last-larval-instar_S1_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/last-larval-instar_S1_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna LLI > data/074_isomiR_reformat/LLI.csv

./074_reformat_isomiRs.pl data/073_isomiR_output/prepupae_S4_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/prepupae_S4_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/prepupae_S4_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/prepupae_S4_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna prepupa > data/074_isomiR_reformat/prepupa.csv 

./074_reformat_isomiRs.pl data/073_isomiR_output/male-pupae_S3_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/male-pupae_S3_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/male-pupae_S3_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/male-pupae_S3_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna male_pupa > data/074_isomiR_reformat/male_pupa.csv

./074_reformat_isomiRs.pl data/073_isomiR_output/female-pupae_S9_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/female-pupae_S9_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/female-pupae_S9_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/female-pupae_S9_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna female_pupa > data/074_isomiR_reformat/female_pupa.csv

./074_reformat_isomiRs.pl data/073_isomiR_output/adult-male_S7_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/adult-male_S7_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/adult-male_S7_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/adult-male_S7_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna male_adult > data/074_isomiR_reformat/male_adult.csv

./074_reformat_isomiRs.pl data/073_isomiR_output/adult-female_S10_L001_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/adult-female_S10_L002_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/adult-female_S10_L003_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna,data/073_isomiR_output/adult-female_S10_L004_R1_001_trim_17-40-N_5_3_adapter.fastq_filterN.mirna female_adult > data/074_isomiR_reformat/female_adult.csv


cat data/074_isomiR_reformat/egg.csv data/074_isomiR_reformat/L1.csv data/074_isomiR_reformat/L5.csv data/074_isomiR_reformat/LLI.csv data/074_isomiR_reformat/prepupa.csv data/074_isomiR_reformat/male_pupa.csv data/074_isomiR_reformat/female_pupa.csv data/074_isomiR_reformat/male_adult.csv data/074_isomiR_reformat/female_adult.csv > data/074_isomiR_reformat/all.csv 
