#!/bin/bash
# script for creation of smallRNA sequences origin distribution

# create output folder
#mkdir -p data/004_condition_smRNA

# build databases
#gunzip -f -k db/tcas5.2_unspliced_transcript.fa.gz
#bwa index db/tcas5.2_unspliced_transcript.fa
#fastq2fasta.pl db/tca_precursor_mirbase_completed_novel.fa > db/tca_precursor_mirbase_completed_novel_dna.fa
#bwa index db/tca_precursor_mirbase_completed_novel_dna.fa
#bwa index db/TCA_all_ncRNA_but_miR.fa_dna.fa
#bwa index db/GCF_000002335.3_Tcas5.2_genomic.fna

# map to genome
for i in data/001_trim_smRNA/* ;
do
	FILEBASENAME=$(basename $i)
	bwa aln -n 1 -o 0 -e 0 -k 1 -t 10 -f data/004_condition_smRNA/${FILEBASENAME}_genome.sai db/GCF_000002335.3_Tcas5.2_genomic.fna ${i};
	bwa samse -f data/004_condition_smRNA/${FILEBASENAME}_genome.sam db/GCF_000002335.3_Tcas5.2_genomic.fna data/004_condition_smRNA/${FILEBASENAME}_genome.sai ${i};
	samtools view -b -f 4 data/004_condition_smRNA/${FILEBASENAME}_genome.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_genome_ual.fq;
	samtools view -b -F 4 data/004_condition_smRNA/${FILEBASENAME}_genome.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_genome_aln.fq;
	#map to ncRNA - map reads that did map to genome
	bwa aln -n 1 -o 0 -e 0 -k 1 -t 10 -f data/004_condition_smRNA/${FILEBASENAME}_ncRNA.sai db/TCA_all_ncRNA_but_miR.fa_dna.fa data/004_condition_smRNA/${FILEBASENAME}_genome_aln.fq;
	bwa samse -f data/004_condition_smRNA/${FILEBASENAME}_ncRNA.sam db/TCA_all_ncRNA_but_miR.fa_dna.fa data/004_condition_smRNA/${FILEBASENAME}_ncRNA.sai data/004_condition_smRNA/${FILEBASENAME}_genome_aln.fq;
	samtools view -b -f 4 data/004_condition_smRNA/${FILEBASENAME}_ncRNA.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_ncRNA_ual.fq;
	samtools view -b -F 4 data/004_condition_smRNA/${FILEBASENAME}_ncRNA.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_ncRNA_aln.fq;
	#map to miRNA - map reads that did not map to ncRNAs
	bwa aln -n 1 -o 0 -e 0 -k 1 -t 10 -f data/004_condition_smRNA/${FILEBASENAME}_miRNA.sai db/tca_precursor_mirbase_completed_novel.fa data/004_condition_smRNA/${FILEBASENAME}_ncRNA_ual.fq;
	bwa samse -f data/004_condition_smRNA/${FILEBASENAME}_miRNA.sam db/tca_precursor_mirbase_completed_novel_dna.fa data/004_condition_smRNA/${FILEBASENAME}_miRNA.sai data/004_condition_smRNA/${FILEBASENAME}_ncRNA_ual.fq;
	samtools view -b -f 4 data/004_condition_smRNA/${FILEBASENAME}_miRNA.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_miRNA_ual.fq;
 	samtools view -b -F 4 data/004_condition_smRNA/${FILEBASENAME}_miRNA.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_miRNA_aln.fq;
	#map to mRNA
	bwa aln -n 1 -o 0 -e 0 -k 1 -t 10 -f data/004_condition_smRNA/${FILEBASENAME}_mRNA.sai db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/${FILEBASENAME}_miRNA_ual.fq;
	bwa samse -f data/004_condition_smRNA/${FILEBASENAME}_mRNA.sam db/tcas5.2_unspliced_transcript.fa data/004_condition_smRNA/${FILEBASENAME}_mRNA.sai data/004_condition_smRNA/${FILEBASENAME}_miRNA_ual.fq;
	samtools view -b -f 4 data/004_condition_smRNA/${FILEBASENAME}_mRNA.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_mRNA_ual.fq;
	samtools view -b -F 4 data/004_condition_smRNA/${FILEBASENAME}_mRNA.sam | bedtools bamtofastq -i - -fq data/004_condition_smRNA/${FILEBASENAME}_mRNA_aln.fq;
	#cleaning up
done
#rm data/004_condition_smRNA/*.sai data/004_condition_smRNA/*.sam 

