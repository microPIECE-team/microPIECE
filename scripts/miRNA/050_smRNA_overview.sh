#!/bin/bash
# script for creation of smallRNA sequences origin distribution

# group reads by condition
mkdir -p data/004_condition_smRNA/egg/
mkdir -p data/004_condition_smRNA/L1/
mkdir -p data/004_condition_smRNA/L5/
mkdir -p data/004_condition_smRNA/LLI/
mkdir -p data/004_condition_smRNA/prepupa/
mkdir -p data/004_condition_smRNA/male_pupa/
mkdir -p data/004_condition_smRNA/female_pupa/
mkdir -p data/004_condition_smRNA/male_adult/
mkdir -p data/004_condition_smRNA/female_adult/

# copy and sort data 
cp data/001_trim_smRNA/egg* data/004_condition_smRNA/egg/
cp data/001_trim_smRNA/First-instar* data/004_condition_smRNA/L1/
cp data/001_trim_smRNA/Fifth-instar* data/004_condition_smRNA/L5/
cp data/001_trim_smRNA/last-larval-instar* data/004_condition_smRNA/LLI/
cp data/001_trim_smRNA/prepupae* data/004_condition_smRNA/prepupa/
cp data/001_trim_smRNA/male-pupae* data/004_condition_smRNA/male_pupa/
cp data/001_trim_smRNA/female-pupae*	data/004_condition_smRNA/female_pupa/
cp data/001_trim_smRNA/adult-male* data/004_condition_smRNA/male_adult/
cp data/001_trim_smRNA/adult-female* data/004_condition_smRNA/female_adult/

# create output folders
mkdir -p data/004_condition_smRNA/egg_quant/
mkdir -p data/004_condition_smRNA/L1_quant/
mkdir -p data/004_condition_smRNA/L5_quant/
mkdir -p data/004_condition_smRNA/LLI_quant/
mkdir -p data/004_condition_smRNA/prepupa_quant/
mkdir -p data/004_condition_smRNA/male_pupa_quant/
mkdir -p data/004_condition_smRNA/female_pupa_quant/
mkdir -p data/004_condition_smRNA/male_adult_quant/
mkdir -p data/004_condition_smRNA/female_adult_quant/

# make databases for alignment
gunzip -f -k db/tcas5.2_unspliced_transcript.fa.gz
bwa index db/tcas5.2_unspliced_transcript.fa
bwa index db/tca_precursor_mirbase_completed_novel.fa
# run quantification
# script 
#	read_folder 
#	genome_BWT_idx 
#	mir_hairpin_BWT_idx 
#	ncRNA_BWT_idx
# 	mRNA_BWT_idx
#	out_folder
#	threads
# 	aligner bwt
#	replicates[4|8]


./051_overview_scan.pl data/004_condition_smRNA/egg/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/egg_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/L1/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/L1_quant/ 100 bwa 8
./051_overview_scan.pl data/004_condition_smRNA/L5/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/L5_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/LLI/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/LLI_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/preupua/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/prepupa_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/male_pupa/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/male_pupa_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/female_pupa/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/female_pupa_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/male_adult/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/male_adult_quant/ 100 bwa 4
./051_overview_scan.pl data/004_condition_smRNA/female_adult/ db/GCF_000002335.3_Tcas5.2_genomic.fna_noWhitespace.fa db/tca_precursor_mirbase_completed_novel.fa db/TCA_all_ncRNA_but_miR.fa_dna.fa db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/female_adult_quant/ 100 bwa 4
