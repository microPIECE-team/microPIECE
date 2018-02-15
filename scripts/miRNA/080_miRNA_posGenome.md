# NAME
080_miRNA_posGenome.pl
# VERSION
0.9
# DEPENDENCIES
- `blast`
- `awk`
# DESCRIPTION
Creates a BLAST database from the species_B genome and runs a blastn search. The results are filtered by 100% identity and coverage hits only. The result contains all genomic loci of all microRNA precursors.
# PARAMETERS
#### makeblastdb
- `-dbtype nucl`
#### blastn
- `-num_threads`
- `-dust no`
- `-soft_masking false`
- `-outfmt 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore`

### awk
- `$3 == 100 && $4 == $5`

# INPUT
- species_B genome.fasta
- species_B_complete_microRNA_precursor_set.fasta

# OUTPUT
blastn output file
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
