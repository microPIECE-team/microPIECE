# NAME
filterBED.sh
# VERSION
0.9
# DEPENDENCIES
`049_bed2signal.pl`
# DESCRIPTION
The shell script calls a perl script that creates `.bed` files. Each `.bed` file represents a different signal strength, from 1, where only one library supports the position, up to the number of libraries being used.
# PARAMETERS
The number of libraries has to be set.
# INPUT
- number of libraries
- `clip_merged.bed` file from previous analysis
# OUTPUT
- `clip_merged_x_of_y_BEDfilter.bed
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
