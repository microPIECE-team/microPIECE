# NAME
blast_qcov_short.pl
# VERSION
0.9
# DEPENDENCIES
- `Getopt::Long`
# DESCRIPTION
Performs a blastn search with included filtering steps
# PARAMETERS
- `-outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq`
- `-word_size 4` := as we search for mature microRNAs, the word size needs to be small
- `-evalue 10000` := short sequences suffer from a high expectation rate of finding the same sequence again by chance
- `-strand plus` := we need only the plus strand for the homology detection
# INPUT
- `-query` := species_B mature microRNAs in multi fasta format 
- `-db` := metazoa matur microRNAs in multi fasta format
- `-out` := output file
- `-num_threads` := number of usable threads
# OUTPUT
blast output in format `6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq`.

The results are filtered in a way that the first ten nucleotides need to be identical. The following nucleotides may contain one mismatch but no gaps.
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
