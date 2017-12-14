#! /bin/bash

time perl 046_merge_bed_files.pl \
     --input 24h_1=../040/SRR5163635_trim_gsnap_piranha_sort.bed \
     --input 24h_2=../040/SRR5163636_trim_gsnap_piranha_sort.bed \
     --input 24h_3=../040/SRR5163637_trim_gsnap_piranha_sort.bed \
     --input 72h_1=../040/SRR5163632_trim_gsnap_piranha_sort.bed \
     --input 72h_2=../040/SRR5163633_trim_gsnap_piranha_sort.bed \
     --input 72h_3=../040/SRR5163634_trim_gsnap_piranha_sort.bed \
     --output clip_merged.bed \
     --log merging_bed_files.log
