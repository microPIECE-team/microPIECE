# NAME
sam2de.pl
# VERSION
0.9
# DEPENDENCIES
- `062_s2d_cfg`
- `Getopt::Long`

# DESCRIPTION
Reads the config file to obtain replicate and condition information:
```
path/to/con1_rep1.sam con1
path/to/con1_rep2.sam con1
path/to/con2_rep1.sam con2
```

# PARAMETERS
none
# INPUT
- `-cfg` := config file `path/to/file.sam condition`
- `-mature_files` := reference mature microRNA set of species_B
# OUTPUT
`.csv` file with `RPM;condition;miRNA`
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
