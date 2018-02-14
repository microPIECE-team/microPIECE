# NAME
quantification.pl
# VERSION
0.9
# DEPENDENCIES
- `miRDeep2`
- `bwa`
- `samtools`
- `061_sam2de.pl`

# DESCRIPTION
# PARAMETERS
#### miRDeep2
- `rna2dna.pl`

#### bwa
- `index`
- `aln`
  - `-n 1`
  - `-o 0`
  - `-e 0`
  - `-k 1`
  - `-t 100`
- `samse`
- `xa2multi.pl`

#### samtools
- `view`
  - `-F 4`


# INPUT
#### rna2dna.pl

#### bwa index 

#### bwa aln

#### bwa samse

#### samtools view

#### xa2multi.pl

#### 061_sam2de.pl

# OUTPUT
#### rna2dna.pl

#### bwa index

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
