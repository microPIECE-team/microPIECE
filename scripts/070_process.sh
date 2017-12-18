#!/bin/bash

REQ_PROGS=(./071_bedtool_discard_sizes.pl bedtools ./072_fasta_uc_and_filter4annotations.pl)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

mkdir -p ../070/

MIN=22
MAX=50

###############################

zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz > ../070/GCF_000004015.4_AaegL3_genomic.fna

# discard sizes
for i in ../050/clip_merged_[0-9]*of6BEDfilter_mapGFF_minLen0.bed
do
    echo "Working on file $i"

    FILEBASENAME=$(basename "$i" .bed)
    ./071_bedtool_discard_sizes.pl "$i" ${MIN} ${MAX} > ../070/${FILEBASENAME}_min${MIN}_max${MAX}.bed

    # sort
    sort -k1,1 -k2,2n ../070/${FILEBASENAME}_min${MIN}_max${MAX}.bed > ../070/${FILEBASENAME}_min${MIN}_max${MAX}_sort.bed

    #get fasta
    bedtools getfasta -s -name -fi \
	     ../070/GCF_000004015.4_AaegL3_genomic.fna \
	     -bed ../070/${FILEBASENAME}_min${MIN}_max${MAX}_sort.bed \
	     -fo ../070/${FILEBASENAME}_min${MIN}_max${MAX}_sort.fa

    # upper case
    ./072_fasta_uc_and_filter4annotations.pl \
	../070/${FILEBASENAME}_min${MIN}_max${MAX}_sort.fa \
	> ../070/${FILEBASENAME}_min${MIN}_max${MAX}_sort_UC.fa
done
