#! /bin/bash
# mapping clip peaks to mRNAs from GFF file
command -v ./1_clip_mapper.pl >/dev/null 2>&1 || { echo "I require 1_clip_mapper.pl but it's not installed. Aborting." >&2; exit 1; }
mkdir ../050/

#copy gff file
zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.gff.gz > ../050/GCF_000004015.4_AaegL3_genomic.gff

# run mapping of peaks to gff mRNAs
./1_clip_mapper.pl ../040/SRR5163632_trim_gsnap_piranha_sort_merge.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/SRR5163632_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
./1_clip_mapper.pl ../040/SRR5163633_trim_gsnap_piranha_sort_merge.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/SRR5163633_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
./1_clip_mapper.pl ../040/SRR5163634_trim_gsnap_piranha_sort_merge.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/SRR5163634_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed

./1_clip_mapper.pl ../040/SRR5163635_trim_gsnap_piranha_sort_merge.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/SRR5163635_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
./1_clip_mapper.pl ../040/SRR5163636_trim_gsnap_piranha_sort_merge.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/SRR5163636_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
./1_clip_mapper.pl ../040/SRR5163637_trim_gsnap_piranha_sort_merge.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
 
