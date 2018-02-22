# NAME
clip2gff.sh
# VERSION
0.9
# DEPENDENCIES
- `zcat`
- `clip_mapper.pl`
# DESCRIPTION
Runs the perl script to assign all transcripts to the CLIP-region for every type of library-support-version.
# PARAMETERS
none
# INPUT
- `<int>` := how many libraries have to support the region
- `species_A.gff` := The GFF of species_A
- `<int>` := minimal length of CLIP-region
- `> out_mapGFF_minLen_int.bed` := output from STDOUT
# OUTPUT
A BED file that has additional information about the transcript IDs in column #4

;annotation=XM_ID1,...,XM_IDn
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
