#!/bin/bash

REQ_PROGS=(./2_bedtool_discard_sizes.pl bedtools ./3_fasta_uc.pl)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done


mkdir ../070/

# discard sizes

./2_bedtool_discard_sizes.pl ../050/clip_merged_4of6BEDfilter_mapGFF_minLen0.bed 22 50 > ../070/clip_merged_mapGFF_minLen0_4of6BEDfilter_min22_max50.bed

# sort
sort -k1,1 -k2,2n ../070/clip_merged_mapGFF_minLen0_4of6BEDfilter_min22_max50.bed > ../070/clip_merged_mapGFF_minLen0_4of6BEDfilter_min22_max50_sort.bed

############################### 

zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz > ../070/GCF_000004015.4_AaegL3_genomic.fna
#get fasta
bedtools getfasta -s -name -fi \
	 ../070/GCF_000004015.4_AaegL3_genomic.fna \
	 -bed ../070/clip_merged_mapGFF_minLen0_4of6BEDfilter_min22_max50_sort.bed \
	 -fo ../070/clip_merged_mapGFF_minLen0_min22_max50_sort.fa

# upper case
./3_fasta_uc.pl \
    ../070/clip_merged_mapGFF_minLen0_4of6BEDfilter_min22_max50_sort.fa \
    > ../070/clip_merged_mapGFF_minLen0_4of6BEDfilter_min22_max50_sort_UC.fa
