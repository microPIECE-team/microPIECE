# NAME
reformat_isomiRs.pl
# VERSION
0.9
# DEPENDENCIES
none
# DESCRIPTION
Takes the `miraligner` output `.mirna` as input. Each condition is run separately. Replicates are provided as comma-separated list.
`miraligner` checks for SNPs, templated additions at 5' and 3' end of the mature microRNA and non-templated additions to the 3' end.
If one type has no entry, the script writes a NULL (SQL-friendly). If all types are NULL, the supporting collapsed read is skipped, as there is obviously no modification.

# PARAMETERS
none
# INPUT
- input list (comma-separated) of replicates of one condition
- condition name
# OUTPUT
microRNA-name;snp;non-template-addition;5prime_addition_or_deletion;3prime_addition_or_deletion;rpm;condition
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
