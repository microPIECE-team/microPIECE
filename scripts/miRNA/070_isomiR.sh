#!/bin/bash
command -v RNAfold >/dev/null 2>&1 || { echo "I require RNAfold but it's not installed. Aborting." >&2; exit 1; }
command -v perl -e 'use RNA::HairpinFigure qw/draw/' >/dev/null 2>&1 || { echo "I require the PERL module RNA::HairpinFigure but it's not installed. Aborting." >&2; exit 1; }
#command -v seqcluster collapse >/dev/null 2>&1 || { echo "I require seqcluster collapse but it's not installed. Aborting." >&2; exit 1; }


mkdir -p data/004_filterN_smRNA/


for i in data/001_trim_smRNA ; do ./071_filter_fastq_N.pl $i > $i.filterN.fastq; done;
mv $i data/004_filterN_smRNA/


wget ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz db/
./072_create_mirbase_struct.pl db/tca_precursor_mirbase_completed_novel.fa db/tca_mature_mirbase_completed_novel.fa
mv custom.str data/070_isomiR_db/
cp db/tca_precursor_mirbase_completed_novel.fa /data/070_isomiR_db/


mkdir -p data/070_isomiR_db/
cp db/ data/070_isomiR_db/
cp .str data/070_isomiR_db/

mkdir -p isomiR_output/

# run miraligner script
./073_seqbuster_pipe.pl data/004_filterN_smRNA/ data/070_isomiR_db/ isomiR_output/ tca


