# NAME
quantification.pl
# VERSION
0.9
# DEPENDENCIES
- `miRDeep2` 
  - `rna2dna.pl`
- `bwa`
  - `xa2multi` 
- `samtools` 
- `061_sam2de.pl` 

# DESCRIPTION
# PARAMETERS
#### bwa
- `index` := indexing the database (the mature microRNAs in our case)
- `aln`
  - `-n 1` := maximal difference
  - `-o 0` := maximum number of gap opens
  - `-e 0` := maximum number of gap extensions
  - `-k 1` := maximum differences in the seed
  - `-t 100` := threads
- `samse` := sai format to sam
- `xa2multi.pl` := belongs to bwa, creates additional lines for multi-mapping reads

#### samtools
- `view` 
  - `-F 4` := only include mapping reads


# INPUT
#### rna2dna.pl
- mature microRNAs of species_B
#### bwa index 
- mature microRNAs of species_B
#### bwa aln
- filtered and trimmed smallRNA data
- mature microRNAs of species_B bwa database
#### bwa samse
- `.sai` output of `bwa aln`
#### samtools view
- `.sam` output of `bwa samse`
#### xa2multi.pl
- mapping reads in `.sam` file of `samtools view`
#### 061_sam2de.pl
- mapping reads with multi-mapping reads, having a separate row in `.sam` file
# OUTPUT
#### rna2dna.pl
- mature microRNAs of species_B in DNA 

#### bwa index
- mature microRNAs of species_B in DNA as BWA database

#### bwa aln
`.sai` format of all reads against the microRNA database

#### bwa samse
`.sam` format of all reads against the microRNA database

#### samtools view
`.sam` format of all aligned reads against the microRNA database

#### xa2multi.pl
`.sam` format of all aligned reads against the microRNA database with every multimapping read in a separate row

#### 061_sam2de.pl
`.csv` format of all microRNAs with RPM expression

First row is header: rpm;condition;miRNA

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
