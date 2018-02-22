# NAME
tarpred.sh
# VERSION
0.9
# DEPENDENCIES
- `zcat`
- `bedtools`
- `095_csv_to_bed.pl`
- `096_mapping.pl`
# DESCRIPTION
Transforms the needle `.csv` output to bedformat, then merges overlapping regions and extracts the transfered CLIP `.fasta` sequence from the transcriptome. Finally a targetprediction on the transfered CLIP-regions with `miranda` is performed.
# PARAMETERS
#### bedtools
- `merge`
- `-c 4 -o collapse`
#### bedtools
- `getfasta`
# INPUT
- set of microRNAs, either from mining-branch of this pipe-line or from other sources
- needle output `.csv` file
- species_B_rna.fa
# OUTPUT
miranda targetprediction file
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
