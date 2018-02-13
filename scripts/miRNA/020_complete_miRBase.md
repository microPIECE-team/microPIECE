# NAME
complete_miRBase.sh
# VERSION
0.9
# DEPENDENCIES
- `021_parse_miRDeep2_output.pl`
# DESCRIPTION
In some cases, miRBase entries are not complete. One arm is e.g. missing. Then sometimes miRBase simply calls the single arm miR-1 instead of miR-1-3p or -5p. This script tries to complete those missing miRbase entries by parsing the mining result of miRDeep2. Another case is that duplicated genomic loci were named e.g. mir-11c-1,mir-11c-2. Their mature sequences on the other hand are named miR-11c-3p/-5p. That means that for two precursor, we only have two instead of four mature sequences. Although they are identical, we need those sequence copies for our analysis. Therefore the user has to provide such cases.

# PARAMETERS
none
# INPUT
- miRdeep2 result in csv format
- `xxx_mature_mirbase.fa`
- genomic precursor copies with identical mature sequences
# OUTPUT
A more complete set of miRNAs in a multi-fasta file.
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
