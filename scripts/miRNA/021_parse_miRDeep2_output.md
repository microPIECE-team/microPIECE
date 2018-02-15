# NAME
parse_miRDeep2_output.pl
# VERSION
0.9
# DEPENDENCIES
- `Getopt::Long`
# DESCRIPTION
It takes the miRDeep2.csv output as source to detect missing annotations in the miRBase multi-fasta file of species_B. The miRBase mature.fa file of species_B is taken as reference. The script can then detect precursors with a missing arm annotation, in case the missing arm was found by miRDeep2. 

The script also takes precursor names with genomic copies in a comma-separated list as input. Those genomic-copy precursors share the same mature sequence and therefore, the script copies the mature sequences and renames them.
Duplicated genomic loci were named e.g. mir-11c-1,mir-11c-2. Their mature sequences on the other hand are named miR-11c-3p/-5p.
This script adds an addtional sequence for 11c-3p and 11c-5p and names the four mature sequences 11c-1-3p, 11c-2-3p, 11c-1-5p and 11-2-5p.
# PARAMETERS
none
# INPUT
- `mirdeep_bwt1.csv` := miRdeep2 output `.csv` file
- `xxx_mature.fa` := mature microRNA multi-fasta fle of species_B
- `<comma_list_of_precursor_copies>` := list of precurosors with genomic copy and identical mature microRNAs
# OUTPUT
A more complete set of mature microRNAs in multi-fasta file.
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
