# NAME
071_filter_fastq_N.pl
# VERSION
0.9
# DEPENDENCIES
none
# DESCRIPTION
If the read contains N's, then it is discarded. Otherwise miraligner would consider N's as SNPs.
# PARAMETERS
none
# INPUT
FASTQ files. 
# OUTPUT
FASTQ files without any N's in their sequence.
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
