# NAME
trimming.sh
# VERSION
0.9
# Dependencies
cutadapt
# DESCRIPTION
Runs cutadapt for each SRA-fastq file and trims the adapter sequence.
# PARAMETERS
`-a GTGTCAGTCACTTCCAGCGG` := adapter
`-m 20` := min. length of retaining reads
`--trim-n` := remove terminal N's from read
`-o out/path/` := path for output files
# INPUT
SRA_raw.fastq files
# OUTPUT
SRA_trimmed.fastq files
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
