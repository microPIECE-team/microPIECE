# NAME
072_create_mirbase_struct.pl
# VERSION
0.9
# DEPENDENCIES
- `RNA::HairpinFigure qw/draw/`
- `Getopt::Long`
- `File::Temp qw(tmpnam)`
- `RNAfold`
# DESCRIPTION
Folds the precursor sequence into a hairpin and format output like miRNA.str.
- Only mature sequences in hairpin are uppercase
- positions of mature sequences in the precursor are printed in the first line
- in case of a missing mature sequence, the miRNA.str entry was taken
# PARAMETERS
#### RNAfold
- `noPS` = no postscript drawing
# INPUT
- `-hairping` := miRNA mature multi-fasta
- `-mature` := miRNA precursor multi-fasta
- `-struct` := structure file from miRBase (miRNA.str)
# OUTPUT
- `custom.str`:= structure file, like miRNA.str, but with the novel microRNAs from previous analysis.
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
