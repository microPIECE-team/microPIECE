# run gsnap to map the reads of CLIP sequencing to AAE genome
command -v gsnap >/dev/null 2>&1 || { echo "I require gsnap but it's not installed. Aborting." >&2; exit 1; }
command -v gmap_build >/dev/null 2>&1 || { echo "I require gsnap but it's not installed. Aborting." >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "I require samtools but it's not installed. Aborting." >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "I require bedtools but it's not installed. Aborting." >&2; exit 1; }



mkdir ../030/
# build genome index for gsnap
zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz > ../030/GCF_000004015.4_AaegL3_genomic.fna
gmap_build -D ../030/gsnap_db -k 15 -d Aedes_aegypti ../030/GCF_000004015.4_AaegL3_genomic.fna


# gsnap mapping
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163635_trim.fastq.gz > ../030/SRR5163635_trim_gsnap.sam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163636_trim.fastq.gz > ../030/SRR5163636_trim_gsnap.sam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163637_trim.fastq.gz > ../030/SRR5163637_trim_gsnap.sam

gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163632_trim.fastq.gz > ../030/SRR5163632_trim_gsnap.sam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163633_trim.fastq.gz > ../030/SRR5163633_trim_gsnap.sam
gsnap --gunzip -N 1 -B 5 --speed 1 -O -A sam -D ../030/gsnap_db -d Aedes_aegypti ../020/SRR5163634_trim.fastq.gz > ../030/SRR5163634_trim_gsnap.sam


# sam2bam 
samtools view -Sb ../030/SRR5163635_trim_gsnap.sam > ../030/SRR5163635_trim_gsnap.bam
samtools view -Sb ../030/SRR5163636_trim_gsnap.sam > ../030/SRR5163636_trim_gsnap.bam
samtools view -Sb ../030/SRR5163637_trim_gsnap.sam > ../030/SRR5163637_trim_gsnap.bam

samtools view -Sb ../030/SRR5163632_trim_gsnap.sam > ../030/SRR5163632_trim_gsnap.bam
samtools view -Sb ../030/SRR5163633_trim_gsnap.sam > ../030/SRR5163633_trim_gsnap.bam
samtools view -Sb ../030/SRR5163634_trim_gsnap.sam > ../030/SRR5163634_trim_gsnap.bam


# bam2bed
bedtools bamtobed -i ../030/SRR5163635_trim_gsnap.bam > ../030/SRR5163635_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163636_trim_gsnap.bam > ../030/SRR5163636_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163637_trim_gsnap.bam > ../030/SRR5163637_trim_gsnap.bed

bedtools bamtobed -i ../030/SRR5163632_trim_gsnap.bam > ../030/SRR5163632_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163633_trim_gsnap.bam > ../030/SRR5163633_trim_gsnap.bed
bedtools bamtobed -i ../030/SRR5163634_trim_gsnap.bam > ../030/SRR5163634_trim_gsnap.bed


