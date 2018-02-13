# NAME
process.sh
# VERSION
0.9
# DEPENDENCIES
- `zcat`
- `bedtool_discard_sizes.pl`
- `sort`
- `bedtools`
- `fasta_uc_and_filter4annotation.pl`
# DESCRIPTION
Discard the `.bed` regions that are smaller than a given minimum and larger than a given maximum. Then extract the nucleotide sequence from the genome of species_A and make all nucleotides uppercase.
# PARAMETERS
#### bedtool_discard_sizes.pl
- `MIN=` := 22 - minimal length of bed region
- `MAX=` := 50 - maximal length of bed region

#### sort
- `-k1,1 -k2,2n`

#### bedtools
- `getfasta`
- `-s`
- `-name`
# INPUT
- `species_A_genome.fa`
- `.bed` files with transcriptome information
# OUTPUT
- length filtered multi-fasta file of CLIP regions.
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
