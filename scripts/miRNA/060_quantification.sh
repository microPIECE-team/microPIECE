#!/bin/bash
mkdir -p data/061_miRNA_expression/
bwa index data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa
# take fq files from data/002_filter_smRNA and map to miRNA mature from data/041_miRDeep_completed_with_novels/
for i in data/002_filter_smRNA/*;
do 
	FILEBASENAME=$(basename $i)
	bwa aln -n 1 -o 0 -e 0 -k 1 -t 100 -f data/061_miRNA_expression/${FILEBASENAME}_mature_miR.sai data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa ${i};
	bwa samse -f data/061_miRNA_expression/${FILEBASENAME}_mature_miR.sam data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa data/061_miRNA_expression/${FILEBASENAME}_mature_miR.sai ${i};
	samtools view -F 4 data/061_miRNA_expression/${FILEBASENAME}_mature_miR.sam -o data/061_miRNA_expression/${FILEBASENAME}_mature_miR_aln.sam;
	./xa2multi.pl data/061_miRNA_expression/${FILEBASENAME}_mature_miR_aln.sam > data/061_miRNA_expression/${FILEBASENAME}_mature_miR_aln_xa2multi.sam;
	rm data/061_miRNA_expression/${FILEBASENAME}_mature_miR.sam data/061_miRNA_expression/${FILEBASENAME}_mature_miR.sai data/061_miRNA_expression/${FILEBASENAME}_mature_miR_aln.sam; 
done

# calculate expression from SAM files 
./061_sam2de.pl -cfg 062_s2d_cfg -mature_file data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa > data/061_miRNA_expression/TCA_miRNA_expression.csv

