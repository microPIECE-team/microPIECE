#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
# result file ends with an empty line

# generate datastructure for expected data 
my $expected_data = qq{rpm;condition;miRNA
250000;a;tca-miR-3811c-1-5p
250000;a;tca-miR-3811c-2-5p
275000;a;tca-miR-1-5p
225000;a;tca-miR-1-3p
200000;b;tca-miR-1-5p
200000;b;tca-miR-1-3p
300000;b;tca-miR-3811c-1-5p
300000;b;tca-miR-3811c-2-5p};



my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/061_sam2de.pl',["-cfg", "060/061_test.cfg", "-mature_file", "060/061_mature.fa"]);
isnt($stdout,"",'061 output is not empty');

#my %got_rpm	= %{&parser($stdout)};
my $got = join("\n",sort(split("\n",$stdout)));
my $exp = join("\n",sort(split("\n",$expected_data)));


is($got,$exp,'061 output as expected');
done_testing();



