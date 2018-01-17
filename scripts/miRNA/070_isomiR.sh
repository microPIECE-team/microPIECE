#!/bin/bash
command -v RNAfold >/dev/null 2>&1 || { echo "I require RNAfold but it's not installed. Aborting." >&2; exit 1; }
command -v java -jar miraligner.jar >/dev/null 2>&1 || { echo "I require miraligner.jar but it's not installed. Aborting." >&2; exit 1; }
command -v perl -e 'use RNA::HairpinFigure qw/draw/' >/dev/null 2>&1 || { echo "I require the PERL module RNA::HairpinFigure but it's not installed. Aborting." >&2; exit 1; }


mkdir -p data/071_filterN_smRNA/
mkdir -p data/072_isomiR_db/

for i in data/001_trim_smRNA/* ; 
do 
	FILEBASENAME=$(basename $i) 
	./071_filter_fastq_N.pl -fq ${i} > data/071_filterN_smRNA/${FILEBASENAME}_filterN.fastq; 
done;


wget -nc ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz -P db/
gunzip -f -k db/miRNA.str.gz
./072_create_mirbase_struct.pl -hairpin data/041_miRDeep_completed_with_novels/tca_precursor_mirbase_completed_novel.fa -mature data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa -struct db/miRNA.str
rm tmp_hairpin.fa tmp_struct.rna

# copy the hairpin and structure file into database folder (necessary for miraligner)
mv custom.str data/072_isomiR_db/miRNA.str
cp data/041_miRDeep_completed_with_novels/tca_precursor_mirbase_completed_novel.fa data/072_isomiR_db/hairpin.fa


mkdir -p data/073_isomiR_output/

# run miraligner script
./073_seqbuster_pipe.pl -fq_path data/071_filterN_smRNA/ -db_path data/072_isomiR_db/ -out_path data/073_isomiR_output/ -species tca


