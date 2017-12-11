#!/bin/bash
REQ_PROGS=(./5_csv_to_bed.pl bedtools ./6_mapping.pl miranda)
for prog in ${REQ_PROGS[*]}
do
    command -v ${prog} >/dev/null 2>&1 || { echo "I require ${prog} but it's not installed. Aborting." >&2; exit 1; }
done



mkdir ../090/ 
zcat ../data/TCA/GCF_000002335.3_Tcas5.2_rna.fna.gz > ../090/GCF_000002335.3_Tcas5.2_rna.fna
cp ../data/TCA/tca_miRNA_mature_dna_novel_add_missing_miRs.fa ../090/tca_miRNA_mature_dna_novel_add_missing_miRs.fa


# csv2bed
./5_csv_to_bed.pl ../080/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle.csv ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle.bed 

#merge transfered clip regions
bedtools merge -i ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle.bed -c 4 -o collapse > ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle_merge.bed

#get fasta sequence of clip regions
bedtools getfasta -name -fi ../090/GCF_000002335.3_Tcas5.2_rna.fna -bed ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle_merge.bed -fo ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle_merge.fa

#targetprediciton
./6_mapping.pl ../090/tca_miRNA_mature_dna_novel_add_missing_miRs.fa ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle_merge.fa ../090/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle_merge_miranda.out
