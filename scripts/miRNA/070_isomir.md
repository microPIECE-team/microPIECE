# NAME
isomir.sh
# VERSION
0.9
# DEPENDENCIES
- `071_filter_fastq_N.pl`
- `072_create_mirbase_struct.pl`
- `073_seqbuster_pipe.pl`
- `074_reformat_isomiRs.pl`
- `wget`
- `gunzip`
- `RNA::HairpinFigure`
- `RNAfold`

# DESCRIPTION
- Filters all reads that contain 'N' characters in the sequence with `071_filter_fastq_N.pl`
- download `miRNA.str.gz` from miRBase
- gunzip `miRNA.str.gz`
- create miRBase structure file, according to `miRNA.str` from miRbase for the novel microRNAs. Use miRNA.str as reference for still incomplete miRNA precursors (missing mature sequence)
- run `miraligner`
- reformat the output to `.CSV` file for each condition and then merge all conditions into one `.CSV` file

# PARAMETERS
none
# INPUT
- trimmed and ncRNA-filtered small RNA sequencing data
- `miRNA.str` microRNA structure file from miRBase
# OUTPUT
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
