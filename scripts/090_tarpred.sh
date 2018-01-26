#!/bin/bash
REQ_PROGS=(./095_csv_to_bed.pl bedtools ./096_mapping.pl miranda)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done

mkdir ../090/
zcat ../data/TCA/GCF_000002335.3_Tcas5.2_rna.fna.gz > ../090/GCF_000002335.3_Tcas5.2_rna.fna
cp miRNA/data/041_miRDeep_completed_with_novels/tca_mature_mirbase_completed_novel.fa ../090/tca_mature_mirbase_completed_novel.fa

for i in ../080/clip_merged_[0-9]*of6BEDfilter_mapGFF_minLen0_min22_max50_sort_UC_needle.csv
do
    echo "Working on file '${i}'"

    FILENAMEBASE=$(basename "${i}" .csv)

    # csv2bed
    ./095_csv_to_bed.pl \
	"${i}" \
	../090/"${FILENAMEBASE}".bed

    #merge transfered clip regions
    bedtools merge -i \
	     ../090/"${FILENAMEBASE}".bed \
	     -c 4 -o collapse \
	     > ../090/"${FILENAMEBASE}"_merge.bed

    #get fasta sequence of clip regions
    bedtools getfasta -name \
	     -fi ../090/GCF_000002335.3_Tcas5.2_rna.fna \
	     -bed ../090/"${FILENAMEBASE}"_merge.bed \
	     -fo ../090/"${FILENAMEBASE}"_merge.fa

    #targetprediciton
    ./096_mapping.pl \
	../090/tca_mature_mirbase_completed_novel.fa \
	../090/"${FILENAMEBASE}"_merge.fa \
	../090/"${FILENAMEBASE}"_merge_miranda.out
done
