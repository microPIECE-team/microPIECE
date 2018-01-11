#!/bin/bash
command -v RNAfold >/dev/null 2>&1 || { echo "I require RNAfold but it's not installed. Aborting." >&2; exit 1; }
command -v perl -e 'use RNA::HairpinFigure qw/draw/' >/dev/null 2>&1 || { echo "I require the PERL module RNA::HairpinFigure but it's not installed. Aborting." >&2; exit 1; }
#command -v seqcluster collapse >/dev/null 2>&1 || { echo "I require seqcluster collapse but it's not installed. Aborting." >&2; exit 1; }


mkdir -p data/004_filterN_smRNA/


for i in data/001_trim_smRNA ; do FILEBASENAME=$(basename $i) ./071_filter_fastq_N.pl data/001_trim_smRNA/{$FILEBASENAME} > data/004_filterN_smRNA/{$FILEBASENAME}_filterN.fastq; done;


wget ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz -P db/
gunzip -f -k db/miRNA.str.gz
./072_create_mirbase_struct.pl db/tca_precursor_mirbase_completed_novel.fa db/tca_mature_mirbase_completed_novel.fa db/miRNA.str
rm tmp_hairpin.fa tmp_struct.rna

# copy the hairpin and structure file into database folder (necessary for miraligner)
mv custom.str data/070_isomiR_db/miRNA.str
cp db/tca_precursor_mirbase_completed_novel.fa /data/070_isomiR_db/hairpin.fa


mkdir -p data/070_isomiR_output/

# run miraligner script
./073_seqbuster_pipe.pl data/004_filterN_smRNA/ data/070_isomiR_db/ data/070_isomiR_output/ tca


