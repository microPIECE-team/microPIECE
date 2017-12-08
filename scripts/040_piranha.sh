#! /bin/bash
command -v Piranha >/dev/null 2>&1 || { echo "I require Piranha but it's not installed. Aborting." >&2; exit 1; }
mkdir ../040/

# performing peak calling with Piranha


Piranha -o ../040/SRR5163632_trim_gsnap_piranha.bed -s ../030/SRR5163632_trim_gsnap.bed
Piranha -o ../040/SRR5163633_trim_gsnap_piranha.bed -s ../030/SRR5163633_trim_gsnap.bed
Piranha -o ../040/SRR5163634_trim_gsnap_piranha.bed -s ../030/SRR5163634_trim_gsnap.bed

Piranha -o ../040/SRR5163635_trim_gsnap_piranha.bed -s ../030/SRR5163635_trim_gsnap.bed
Piranha -o ../040/SRR5163636_trim_gsnap_piranha.bed -s ../030/SRR5163636_trim_gsnap.bed
Piranha -o ../040/SRR5163637_trim_gsnap_piranha.bed -s ../030/SRR5163637_trim_gsnap.bed



# sorting
sort -k1,1 -k2,2n ../040/SRR5163632_trim_gsnap_piranha.bed > ../040/SRR5163632_trim_gsnap_piranha_sort.bed
sort -k1,1 -k2,2n ../040/SRR5163633_trim_gsnap_piranha.bed > ../040/SRR5163633_trim_gsnap_piranha_sort.bed 
sort -k1,1 -k2,2n ../040/SRR5163634_trim_gsnap_piranha.bed > ../040/SRR5163634_trim_gsnap_piranha_sort.bed

sort -k1,1 -k2,2n ../040/SRR5163635_trim_gsnap_piranha.bed > ../040/SRR5163635_trim_gsnap_piranha_sort.bed
sort -k1,1 -k2,2n ../040/SRR5163636_trim_gsnap_piranha.bed > ../040/SRR5163636_trim_gsnap_piranha_sort.bed
sort -k1,1 -k2,2n ../040/SRR5163637_trim_gsnap_piranha.bed > ../040/SRR5163637_trim_gsnap_piranha_sort.bed

# merge
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163632_trim_gsnap_piranha_sort.bed > ../040/SRR5163632_trim_gsnap_piranha_sort_merge.bed
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163633_trim_gsnap_piranha_sort.bed > ../040/SRR5163633_trim_gsnap_piranha_sort_merge.bed 
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163634_trim_gsnap_piranha_sort.bed > ../040/SRR5163634_trim_gsnap_piranha_sort_merge.bed

bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163635_trim_gsnap_piranha_sort.bed > ../040/SRR5163635_trim_gsnap_piranha_sort_merge.bed
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163636_trim_gsnap_piranha_sort.bed > ../040/SRR5163636_trim_gsnap_piranha_sort_merge.bed
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163637_trim_gsnap_piranha_sort.bed > ../040/SRR5163637_trim_gsnap_piranha_sort_merge.bed
