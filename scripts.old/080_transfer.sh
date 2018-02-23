#!/bin/bash
REQ_PROGS=(./081_map_clip_gff_needle.pl ./085_parse_gff_return_longest_transcript.pl)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

mkdir -p ../080/

# filter GFF files for unqiue longest transcripts
zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.gff.gz > ../080/GCF_000004015.4_AaegL3_genomic.gff
zcat ../data/TCA/GCF_000002335.3_Tcas5.2_genomic.gff.gz > ../080/GCF_000002335.3_Tcas5.2_genomic.gff
# AAE
./085_parse_gff_return_longest_transcript.pl ../080/GCF_000004015.4_AaegL3_genomic.gff > ../080/GCF_000004015.4_AaegL3_genomic_XM_XP_unique.csv 2> ../080/GCF_000004015.4_AaegL3_genomic_XM_XP_unique.err
# TCA
./085_parse_gff_return_longest_transcript.pl ../080/GCF_000002335.3_Tcas5.2_genomic.gff > ../080/GCF_000002335.3_Tcas5.2_genomic_XM_XP_unique.csv 2> ../080/GCF_000002335.3_Tcas5.2_genomic_XM_XP_unique.err

zcat ../data/TCA/GCF_000002335.3_Tcas5.2_rna.fna.gz > ../080/GCF_000002335.3_Tcas5.2_rna.fna

for i in ../070/clip_merged_[0-9]*of6BEDfilter_mapGFF_minLen0_min22_max50_sort_UC.fa
do
    echo "Working on file ${i}"

    FILEBASENAME=$(basename "${i}" .fa)

    # transfer with needle
    ./081_map_clip_gff_needle.pl \
	../080/GCF_000004015.4_AaegL3_genomic_XM_XP_unique.csv \
	../010/TCA_vs_AAE.proteinortho \
	${i} \
	../080/GCF_000002335.3_Tcas5.2_genomic_XM_XP_unique.csv \
	../080/GCF_000002335.3_Tcas5.2_rna.fna \
	../080/${FILEBASENAME}_needle.csv \
	> ../080/${FILEBASENAME}_needle.aln
done
