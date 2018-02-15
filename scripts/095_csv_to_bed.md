# NAME
csv_to_bed.pl
# VERSION
0.9
# DEPENDENCIES
- `File::Temp qw(tmpnam)`
# DESCRIPTION
Parses the needle `.csv` output file and transforms it into `.bed` format for further usage.
# PARAMETERS
none
# INPUT
- needle `.csv` output file
# OUTPUT
- needle `.csv` output file in `.bed` format

The header section of the `.bed` file contains information about the previous analysis, earlier explained.
  - `QRY`
  - `TAR`
  - `CLIP`
  - `IDENTITY`
  - `COVERAGE`
  - `QRY_GAPS`
  - `TAR_GAPS`
  - `TOTAL_GAPS`
  - `MM`
  - `MATCHES`
  - `TAR_START`
  - `TAR_STOP`
 
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
