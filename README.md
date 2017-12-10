## tca_miRNA_data_generation
All stuff which is required to generate data for Daniels miRNA project

# 005_getFiles.sh
Download all files needed for AGO CLIP transfer of AAE to TCA.
Needs fastq-dump from the NCBI SRA Toolkit and wget

# 010_proteinortho.sh
Computes the orthologous proteins of AAE and TCA
Needs ProteinOrtho and BLAST+

# 020_trimming.sh
Trimming of the AGO CLIP sequencing data
Needs cutadapt

# 025_build_db.sh
Builds the GSNAP database out of the AAE genome
Needs gmap-gsnap

# 030_clip_mapping.sh
Mapping of the AGO CIP reads against the AAE genome
Needs gmap-gsnap

# 040_piranha.sh
Peak Calling of AGO binding regions
Needs Piranha

# 045_bedtools_merge.sh
Merges Peaks, result should stay the same
Needs Bedtools

# 050_clip2gff.sh
Maps the Peak regions to the mRNAs of AAE

# 060_intersect.sh
Merges the replicates

# 070_process.sh
Merges the two conditions

# 080_transfer.sh
Extraction of the longest transcripts for each gene from the GFF files of AAE and TCA
Transfer of the AAE CLIP regions to the orthologous longest transcripts from TCA


