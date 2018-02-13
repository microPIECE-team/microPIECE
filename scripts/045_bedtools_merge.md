# NAME
bedtools_merge.sh
# VERSION
0.9
# DEPENDENCIES
`046_merge_bed_files.pl`
# DESCRIPTION
Calls the perl script to merge the CLIP signal `.bed` files into one. The `.bed` file contains a custom information in column #4. There it shows how many libraries support the genomic location as AGO binding region.
# PARAMETERS
- `--input conX_repY=/path/to/SRR_trim_gsnap_piranha_sort.bed` := for each library 
# INPUT
`/path/to/SRR_trim_gsnap_piranha_sort.bed`
# OUTPUT
`clip_merged.bed`
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
