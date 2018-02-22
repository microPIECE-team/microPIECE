# NAME
seq_buster_pipe.pl
# VERSION
0.9
# DEPENDENCIES
- `miraligner`
- `java`
# DESCRIPTION
Collapse `fastq` file and run `miraligner`
# PARAMETERS
#### miraligner
- `-sub 1` := mismatches allowed
- `-trim 3` := max allowed nt for trimming 
- `-add 3` := max allowed nt for addition
- `-freq` := add frequency to the output
- `-s` := species 3 letter code
# INPUT
trimmed and filtered `.fastq` files
# OUTPUT
`.mirna` file with results from miraligner
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
