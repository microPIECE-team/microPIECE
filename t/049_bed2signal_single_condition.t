#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

use bed;

my $first_file = qq{NC_010241.1	0	1	length=1;counts=1/1	.	+
NC_010241.1	30	45	length=15;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,2/2,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	10	20	length=10;counts=2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2	.	+
NC_010241.1	10	20	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	-
NC_010241.1	150	165	length=15;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	70	80	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+};

my $expected_output = qq{NC_010241.1	10	20	length=10;counts=2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2;originallength=10;originalbedline=3;subregion=0,9	.	+
NC_010241.1	39	40	length=1;counts=2/2;originallength=15;originalbedline=2;subregion=9,9	.	+
};

my ($fh, $filename) = File::Temp::tempfile();

print $fh $first_file;
close($fh) || die;

my ($return,$stdout,$stderr)=run_script('../scripts/049_bed2signal.pl');
is(Test::Script::Run::last_script_exit_code(), 2, 'Without arguments is should exit with exit code 2');

($return,$stdout,$stderr)=run_script('../scripts/049_bed2signal.pl' ,[$filename, 2]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With single file/condition is should exit with exit code 0');

my $got	= &bed::parser($stdout);
my $expected = &bed::parser($expected_output);
is_deeply($got,$expected,'output as expected');

done_testing();
