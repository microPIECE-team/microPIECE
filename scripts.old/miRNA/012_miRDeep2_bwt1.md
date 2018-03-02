# NAME
miRDeep2_bwt1.pl
# VERSION
0.9
# DEPENDENCIES
- `miRDeep2.pl`
  - `fastq2fastq.pl`
  - `remove_white_space_in_id.pl`
  - `collapse_reads_md.pl`
  - `mapper.pl`
- `Getopt::Long`
- bowtie
# DESCRIPTION
Translates the miRNA sequences from RNA to DNA, removes whitespaces in the headers, collapses the reads 
# PARAMETERS
#### bowtie
- `build` := creates a bowtie1 database out of the genome

#### mapper.pl
- `-c` := input is in fasta format
- `-q` := map with one mismatch in seed
- `-n` := overwrite existing files
- `-l 17` := discard reads shorter than 17
- `-p` := map to genome bowtie index

#### miRDeep2.pl
- `-P` := use 3p-5p notation

# INPUT
#### miRDeep2.pl
- `reads` := fasta reads
- `genome` := multifasta file of genome, unique identifiers needed
- `mapped reads in arf format` := taken from mapper.pl
- `reference_miRNAs` := species mature miRNAs
- `other_miRNAs` := other mature microRNAs than from species
- `reference_precursors` := precursors of species
# OUTPUT
miRDeep2 output in CSV format.
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
