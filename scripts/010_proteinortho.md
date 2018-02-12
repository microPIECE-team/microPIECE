# NAME
010_proteinortho.pl
# VERSION
0.9
# Dependencies
- makeblast
- proteinortho5.pl
# DESCRIPTION
Computing blastp databases of both protein.fa files from species A and B. Then proteinortho detects orthologous genes between those two species.
Proteinortho automatically uses all available CPUs, but can be limited with `-cpus=` .
# PARAMETERS
`-clean`

`--project=output/folder/`
# INPUT
- species_A_protein.fa
- species_B_protein.fa
# OUTPUT
Tab-separated .proteinortho file:
- `Species` := number of species in this orthology group
- `Genes` := number of genes in this orthology group. Genes>Species indicates co-orthologs.
- `Alg.-Conn.` := algebraic connectivity
- `species_A_protein.fa` := proteinIDs involved from species_A
- `species_B_protein.fa` := proteinIDs involved from species_B
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
