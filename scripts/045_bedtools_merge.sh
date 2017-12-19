#! /bin/bash
# check for required programs
command -v perl -e 'Pod::Usage' >/dev/null 2>&1 || { echo "I require the perl module Pod::Usage but it's not installed. Aborting." >&2; exit 1; }
command -v perl -e 'Log::Log4perl qw(:easy)' >/dev/null 2>&1 || { echo "I require the perl module Log::Log4perl qw(:easy) but it's not installed. Aborting." >&2; exit 1; }
command -v perl -e 'Getopt::Long' >/dev/null 2>&1 || { echo "I require the perl module Getopt::Long but it's not installed. Aborting." >&2; exit 1; }
command -v perl -e 'Data::Dumper;' >/dev/null 2>&1 || { echo "I require the perl module Data::Dumper but it's not installed. Aborting." >&2; exit 1; }



time perl 046_merge_bed_files.pl \
     --input 24h_1=../040/SRR5163635_trim_gsnap_piranha_sort.bed \
     --input 24h_2=../040/SRR5163636_trim_gsnap_piranha_sort.bed \
     --input 24h_3=../040/SRR5163637_trim_gsnap_piranha_sort.bed \
     --input 72h_1=../040/SRR5163632_trim_gsnap_piranha_sort.bed \
     --input 72h_2=../040/SRR5163633_trim_gsnap_piranha_sort.bed \
     --input 72h_3=../040/SRR5163634_trim_gsnap_piranha_sort.bed \
     --output ../040/clip_merged.bed \
     --log ../040/merging_bed_files.log
