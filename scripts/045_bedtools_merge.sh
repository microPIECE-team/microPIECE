#! /bin/bash
# check for required programs
for i in Data::Dumper Pod::Usage Log::Log4perl Getopt::Long
do
perl -M${i} -e 'print "YEAH\n"' &>/dev/null;
if [ $? = 0 ]; then echo "Module ${i} installed!"; else echo "Module ${i} missing!"; exit 1; fi
done



#time perl 046_merge_bed_files.pl \
#     --input 24h_1=../040/SRR5163635_trim_gsnap_piranha_sort.bed \
#     --input 24h_2=../040/SRR5163636_trim_gsnap_piranha_sort.bed \
#     --input 24h_3=../040/SRR5163637_trim_gsnap_piranha_sort.bed \
#     --input 72h_1=../040/SRR5163632_trim_gsnap_piranha_sort.bed \
#     --input 72h_2=../040/SRR5163633_trim_gsnap_piranha_sort.bed \
#     --input 72h_3=../040/SRR5163634_trim_gsnap_piranha_sort.bed \
#     --output ../040/clip_merged.bed \
#     --log ../040/merging_bed_files.log
