# NAME
bedtool_discard_sizes.pl
# VERSION
0.9
# DEPENDENCIES
none
# DESCRIPTION
Simply checks for the length of the region and discards it, if it is longer or shorter than the given thresholds.
# PARAMETERS
# INPUT
- `.bed` := bed file of previous analysis
- `<int>` := minimal length
- `<int>` := maximal length
# OUTPUT
The bed file that only contains regions between minimal and maximal lengths.
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
