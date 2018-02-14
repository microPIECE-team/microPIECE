# NAME
curated_mirdeep2fasta.pl
# VERSION
0.9
# DEPENDENCIES
- `Getopt::Long`
# DESCRIPTION
Extracts the novel microRNA precursors and mature sequences from miRDeep2 output. Replaces U with T and makes all nucleotides uppercase.
# PARAMETERS
none
# INPUT
- `-csv` := miRDeep2 `.csv` output file
- `-cutoff` := minimal score in miRDeep2, to be considered as novel microRNA 
# OUTPUT
Complete precursor and mature sets of novel and existing microRNAs. 
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
