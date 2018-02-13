# NAME
clip_mapping.sh
# VERSION
0.9
# DEPENDENCIES
- `gsnap`
- `samtools`
- `bedtools`
# DESCRIPTION
# PARAMETERS
`gsnap` :
- `--gunzip`
- `-N 1`
- `-B 5`
- `--speed 1`
- `-O`
- `-A sam`
- `-t`
- `-D`
- `-d`
`samtools`:
- `view -Sb`
- `sort -o`
`bedtools`:
- `bamtobed`
# INPUT
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
