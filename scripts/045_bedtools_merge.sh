#! /bin/bash
# merge
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163632_trim_gsnap_piranha_sort.bed > ../040/SRR5163632_trim_gsnap_piranha_sort_merge.bed
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163633_trim_gsnap_piranha_sort.bed > ../040/SRR5163633_trim_gsnap_piranha_sort_merge.bed 
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163634_trim_gsnap_piranha_sort.bed > ../040/SRR5163634_trim_gsnap_piranha_sort_merge.bed

bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163635_trim_gsnap_piranha_sort.bed > ../040/SRR5163635_trim_gsnap_piranha_sort_merge.bed
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163636_trim_gsnap_piranha_sort.bed > ../040/SRR5163636_trim_gsnap_piranha_sort_merge.bed
bedtools merge -s -c 4,5,6 -o distinct -i ../040/SRR5163637_trim_gsnap_piranha_sort.bed > ../040/SRR5163637_trim_gsnap_piranha_sort_merge.bed
