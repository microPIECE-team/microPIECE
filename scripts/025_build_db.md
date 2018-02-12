# NAME
build_db.sh
# VERSION
0.9
# Dependencies
- gmap_build
# DESCRIPTION
Build an alignment database for `gsnap` from the species_A genome.
# PARAMETERS
- `-k 15` := kmer size 15 (requires 4GB of RAM)
- `-g` := input is gzipped
- `-D` := Destination folder
- `-d Species_name` := Genome name
# INPUT
- `species_A_genome.fa`
- `Species_name`
# OUTPUT
Indexed genome of species_A, ready to use for gsnap alignment.
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
