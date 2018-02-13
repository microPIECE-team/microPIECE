# NAME
bed2signal.pl
# VERSION
0.9
# DEPENDENCIES
none
# DESCRIPTION
Takes the `.bed` with the custom #4 column as input, together with the wanted signal strength (how many libraries have to support the genomic loci).
# PARAMETERS
none
# INPUT
- `file.bed`
- `<int_signal_strength>`
# OUTPUT
A `.bed` file that only contains those regions that have at least the given number of supportive libraries.

Adds additional information to column #4 : 

`length=<int>;subregion=<int>,<int>;originallength=<int>;originalbedline=<int>`
  - length := new length of region
  - subregion := string-substring indicator
  - originallength := original length of bed region from where this subregion was taken
  - originalbedline := the line number of the full-length region in the original bed file
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
