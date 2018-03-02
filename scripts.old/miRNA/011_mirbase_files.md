# NAME
mirbase_files.pl
# VERSION
0.9
# DEPENDENCIES
- `Getopt::Long`
# DESCRIPTION
Takes the `mature.fa`, `hairpin.fa` and `organism.txt` files from miRBase.org together with the 3letter code of the species and divides the files into mature sequences of the species, precursor sequences of the species and all other mature sequences. MicroRNAs other than metaoza are filtered out.
# PARAMETERS
none
# INPUT
- `organism` := organism.txt
- `species` := 3letter code of species of interest (species_B in our case)
- `precursor_file` := hairpin.fa
- `mature_file` := mature.fa
# OUTPUT
- `xxx_mature_sequences` := all mature miRNAs from the species 
- `all_other_mature_sequences` := all other miRNAs (metazoa only)
- `xxx_precursor_sequences` := all precursor miRNAs from the species
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
