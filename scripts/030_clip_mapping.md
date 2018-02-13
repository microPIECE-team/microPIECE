# NAME
clip_mapping.sh
# VERSION
0.9
# DEPENDENCIES
- `gsnap`
- `samtools`
- `bedtools`
# DESCRIPTION
# PARAMETERS
`gsnap` :
- `--gunzip` := input needs to be unzipped
- `-N 1` := look for novel splice sites
- `-B 5` := batch mode 5
- `--speed 1` := slowest, but highest accuarcy
- `-O` := Output is in the same order as input
- `-A sam` := change output type to sam
- `-t` := threads to be used
- `-D` := genome directory
- `-d` := genome database

`samtools`:
- `view -Sb` := change sam to bam 
- `sort -o` := sort bam file

`bedtools`:
- `bamtobed` := change bam to bed
# INPUT
- `SRR_trim.fastq.gz`
- `species_A database`
# OUTPUT
- `SRR_trim_gsnap.bed`
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
