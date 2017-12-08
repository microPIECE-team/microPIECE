#! /bin/bash
# merge replicates by intersecting regions
command -v bedtools >/dev/null 2>&1 || { echo "I require bedtools but it's not installed. Aborting." >&2; exit 1; }

mkdir ../060/

bedtools intersect -s -a ../050/SRR5163632_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed -b ../050/SRR5163633_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed ../050/SRR5163634_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed > ../060/SRR5163632_SRR5163633_SRR5163634_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
 
bedtools intersect -s -a ../050/SRR5163635_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed -b ../050/SRR5163636_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed ../050/SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed > ../060/SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0.bed
