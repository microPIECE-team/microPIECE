#!/bin/bash

REQ_PROGS=(./2_bedtool_discard_sizes.pl bedtools ./3_fasta_uc.pl)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done


mkdir ../070/

# discard sizes

./2_bedtool_discard_sizes.pl ../060/SRR5163632_SRR5163633_SRR5163634_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect.bed 22 50 > ../070/SRR5163632_SRR5163633_SRR5163634_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50.bed
./2_bedtool_discard_sizes.pl ../060/SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect.bed 22 50 > ../070/SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50.bed


# cat 
cat ../070/SRR5163632_SRR5163633_SRR5163634_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50.bed ../070/SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50.bed > ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat.bed

# sort
sort -k1,1 -k2,2n ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat.bed > ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort.bed

# merge
bedtools merge -s -c 4,5,6 -o distinct -i ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort.bed > ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge.bed


############################### 

zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz > ../070/GCF_000004015.4_AaegL3_genomic.fna
#get fasta
bedtools getfasta -s -name -fi ../070/GCF_000004015.4_AaegL3_genomic.fna -bed ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge.bed -fo ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge.fa

# upper case
./3_fasta_uc.pl ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge.fa > ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC.fa

