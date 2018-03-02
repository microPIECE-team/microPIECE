# NAME
map_clip_gff_needle.pl
# VERSION
0.9
# DEPENDENCIES
- `File::Temp qw(tmpnam)`
- `needle`

# DESCRIPTION
The CLIP fasta sequences from species_A are aligned with needle against the orthologous transcripts from species_B.
# PARAMETERS
#### needle
- `-datafile EDNACUSTOM`
- `-endweight Y`
- `-gapopen 5`
- `-gapextend 2`
- `-auto`
- `-aformat markx3`

# INPUT
- `species_A_gff.csv` := transcript_ID and protein_ID csv file of species_A
- `species_A_vs_species_B.proteinortho.pl` := proteinortho result file of species_A vs species_B
- `clip.bed` := CLIP fasta file
- `species_B_gff.csv` := transcript_ID and protein_ID csv file of species_B
- `outfile` := output file
# OUTPUT
A `.csv` file with the first line, being the header (without #). Sections are tab-separated.
- `XM_QRY;XP_QRY` := species_A_transcript_ID;species_A_protein_ID
||      
- `XM_TAR;XP_TAR` := species_B_transcript_ID;species_A_protein_ID
||     
- `qry_header tar_header identity coverage QRY_len TAR_len qry_gaps tar_gaps total_gaps MM matches start stop seq` :=
  - header of CLIP sequence
  - species_B_transcriptome header
  - needle identity
  - needle coverage
  - needle species_A clip sequence length
  - needle species_B transcript sequence length
  - needle species_A clip sequence gaps
  - needle species_B transcript sequence gaps
  - needle sum of clip and transcript sequence gaps
  - needle mismatches in alignment
  - needle matches in alignment
  - needle alignment start in transcript
  - needle alignment stop in transcript
  - needle alignment sequence
  
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
