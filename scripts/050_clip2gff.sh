#! /bin/bash
# mapping clip peaks to mRNAs from GFF file
command -v ./051_clip_mapper.pl >/dev/null 2>&1 || { echo "I require 1_clip_mapper.pl but it's not installed. Aborting." >&2; exit 1; }
mkdir -p ../050/

#copy gff file
zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.gff.gz > ../050/GCF_000004015.4_AaegL3_genomic.gff

# run mapping of peaks to gff mRNAs
for i in ../040/clip_merged_*of6BEDfilter.bed;
do
    echo "Working on file $i"

    ./051_clip_mapper.pl "$i" ../050/GCF_000004015.4_AaegL3_genomic.gff 0 >../050/$(basename "$i" .bed)_mapGFF_minLen0.bed
done
