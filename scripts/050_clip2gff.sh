#! /bin/bash
# mapping clip peaks to mRNAs from GFF file
command -v ./1_clip_mapper.pl >/dev/null 2>&1 || { echo "I require 1_clip_mapper.pl but it's not installed. Aborting." >&2; exit 1; }
mkdir ../050/

#copy gff file
zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.gff.gz > ../050/GCF_000004015.4_AaegL3_genomic.gff

# run mapping of peaks to gff mRNAs
./1_clip_mapper.pl ../040/clip_merged.bed ../050/GCF_000004015.4_AaegL3_genomic.gff 0 > ../050/clip_merged_mapGFF_minLen0.bed
 
