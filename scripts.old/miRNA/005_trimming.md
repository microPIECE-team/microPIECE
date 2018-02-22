# NAME
trimming.sh
# VERSION
0.9
# DEPENDENCIES
- `gunzip`
- `bwa`
- `samtools`
# DESCRIPTION
Uses `cutadapt` to trim the raw smallRNA seqeuncing reads and filters them against a set of ncRNAs (miRNAs excluded). 
# PARAMETERS
#### bwa
- `index`

- `aln`
  - `n 1` := maximal difference
  - `o 0` := maximum number of gap opens
  - `e 0` := maximum number of gap extensions
  - `k 1` := maximum differences in the seed
  - `t 100` := threads
  
- `samse`

#### samtools
- `view`
  - `-b` := output in bam format
  - `-f 4` := only include unaligned reads
# INPUT
#### cutadapt
- adapter sequence
- smallRNA raw reads

#### bwa 
- species_B ncRNA set

# OUTPUT
Trimmed and ncRNA-filtered smallRNA sequencing reads.

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
