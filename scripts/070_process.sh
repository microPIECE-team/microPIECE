#!/bin/bash

mkdir ../070/

# discard sizes

./2_bedtool_discard_sizes.pl ../060/.bed 50 > ../070/_max50.bed


# cat 

# sort
sort -k1,1 -k2,2n

# merge
bedtools merge -s -c 4,5,6 -o distinct -i


############################### 


#get fasta
bedtools getfasta -s -name -fi -bed -fo

# upper case
./3_fasta_uc.pl

