# NAME
mining.sh
# VERSION
0.9
# DEPENDENCIES
- `gunzip`
- `011_mirbase_files.pl`
- `012_miRDeep2_bwt1.pl`
# DESCRIPTION
Downloads the current miRBase files from hairpin, mature and organism.

# PARAMETERS
none
# INPUT
#### `011_mirbase_files.pl`
- `-species` : = 3letter code of species of interest (species_B in our case)
- `-precursor_file` := hairpin.fa
- `-mature_file` := mature.fa
- `-organism` := organism.txt

#### `012_miRDeep2_bwt1.pl`
- `-dir` := input directory
- `-out` := output directory
- `-ref_genome` := species_A_genome.fa
- `-species_mature_miRs` := xxx_mature_mirbase.fa
- `-other_mature_miRs` := mature.fa-no-xxx.fa
- `-species_precursor_mirs` := xxx_precursor_mirbase.fa
- `threads` := number of threads to use

# OUTPUT
#### `011_mirbase_files.pl`
- `xxx_mature_sequences` := all mature miRNAs from the species 
- `all_other_mature_sequences` := all other miRNAs (metazoa only)
- `xxx_precursor_sequences` := all precursor miRNAs from the species

#### `012_miRDeep2_bwt1.pl`
- `result_bwt1.csv` := miRDeep2 output file

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
