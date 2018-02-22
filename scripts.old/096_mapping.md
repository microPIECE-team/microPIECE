# NAME
mapping.pl
# VERSION
0.9
# DEPENDENCIES
- `File::Temp qw(tmpnam)`
- `miranda`
# DESCRIPTION
Uses the miRNA set, mined in the other branch of this pipeline or a foreing set of microRNAs to predict their targets on the transfered CLIP-regions.
# PARAMETERS
#### miranda
- `quiet`
# INPUT
- `<mature_miRNAs.fa>` := mined microRNAs from other branch of the pipeline or foreign set of mature microRNAs
- `<transfered_clip_sequences.fa>` := A multi-fasta file of transfered CLIP sequences 
- `<outfile>` := output file
# OUTPUT
Tab separated result file with header line

- query_miRNA := microRNA name
- target_mRNA := header of transfered CLIP region (very long, as there is all information from needle inside)
- score := miranda scoring
- kcal/mol := binding energy
- Query-aln_start := start of alignment at query sequence
- Query-aln_end := stop of alignment at query sequence
- Subject-aln_start := start of alignment at target sequence
- Subject_aln-End := stop of alignment at target sequence
- Al-Len := length of alignment
- Subject-Identity := sequence identity in CLIP sequence
- Query-Identity := sequence identity in microRNA sequence

# CHANGELOG
- 2018-02-12 Release version 0.9
# KNOWN BUGS
no known bugs
# LICENSE
This program is released under GPLv2. For further license information, see LICENSE.md shipped with this program.
Copyright(c)2018 Daniel Amsel and Frank Förster (employees of Fraunhofer Institute for Molecular Biology and Applied Ecology IME).
All rights reserved.
# CONTACT
daniel.amsel@ime.fraunhofer.de
frank.foerster@ime.fraunhofer.de
