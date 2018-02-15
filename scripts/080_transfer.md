# NAME
transfer.sh
# VERSION
0.9
# DEPENDENCIES
- `zcat`
- `parse_gff_return_longest_transcript.pl`
- `map_clip_gff_needle.pl`
# DESCRIPTION
Calls the GFF parser `parse_gff_return_longest_transcript.pl`, to obtain the ID longest transcript and its protein ID as CSV.

Then `map_clip_gff_needle.pl` uses `needle` and `EDNACUSTOM` to transfer the CLIP regions to orthologous transcriptomes by a needle-alignment.
# PARAMETERS
none
# INPUT
- `species_A.gff` := `.gff` file of speices_A
- `species_B.gff` := `.gff` file of species_B
- `species_A_vs_species_B.proteinortho` := proteinortho result file
- `species_B_rna.fa` := transcriptome of species_B
- `EDNACUSTOM`  := NEEDLE scoring file, modified
# OUTPUT
- `needle_output.csv`
- `needle_alignment.aln`
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
