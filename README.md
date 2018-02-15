# tca_miRNA_data_generation
All stuff which is required to generate data for Daniels miRNA project

## How to run
### First, we run the microRNA pipeline
`mkdir -p scripts/miRNA/data/000_raw_smRNA/`<br />
Copy smallRNA sequencing raw fastq files into the folder and run the pipeline with <br />
`scripts/miRNA/run_all.sh`<br />
### Second, we do the RNA-seq analysis

### Third, we perform the CLIP transfer


# Script details
### 005_getFiles.sh
Download all files needed for AGO CLIP transfer of AAE to TCA.
Needs fastq-dump from the NCBI SRA Toolkit
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

### 010_proteinortho.sh
Computes the orthologous proteins of AAE and TCA
Needs ProteinOrtho and BLAST+
https://www.bioinf.uni-leipzig.de/Software/proteinortho/
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

### 020_trimming.sh
Trimming of the AGO CLIP sequencing data
Needs cutadapt
https://github.com/marcelm/cutadapt

### 025_build_db.sh
Builds the GSNAP database out of the AAE genome
Needs gmap-gsnap
http://research-pub.gene.com/gmap/

### 030_clip_mapping.sh
Mapping of the AGO CIP reads against the AAE genome
Needs gmap-gsnap, samtools and bedtools
http://research-pub.gene.com/gmap/
http://samtools.sourceforge.net/
http://bedtools.readthedocs.io/en/latest/

### 040_piranha.sh
Peak Calling of AGO binding regions
Needs Piranha
http://smithlabresearch.org/software/piranha/

### 045_bedtools_merge.sh
Merges Peaks, result should stay the same
Needs Bedtools
http://bedtools.readthedocs.io/en/latest/

### 050_clip2gff.sh
Maps the Peak regions to the mRNAs of AAE

### 060_intersect.sh
Merges the replicates
Needs Bedtools
http://bedtools.readthedocs.io/en/latest/

### 070_process.sh
Merges the two conditions
Needs Bedtools
http://bedtools.readthedocs.io/en/latest/

### 080_transfer.sh
Extraction of the longest transcripts for each gene from the GFF files of AAE and TCA
Transfer of the AAE CLIP regions to the orthologous longest transcripts from TCA
Needs needle from the EMBOSS package
ftp://emboss.open-bio.org/pub/EMBOSS/


### 090_tarpred.sh
Extracts the transfered CLIP regions in fasta format and predicts the miRNA targets
Needs bedtools and miranda
http://bedtools.readthedocs.io/en/latest/
http://34.236.212.39/microrna/getDownloads.do

