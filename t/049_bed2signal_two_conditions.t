#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

use bed;

my $first_file = qq{NC_010241.1	10	20	length=10;counts=2/1/1,1/1/0,3/1/2,1/1/0,5/2/3,1/0/1,10/8/2,1/1,8/1/7,1/1	.	-
};

my $expected_output = qq{NC_010241.1	10	11	length=1;counts=2/1/1;originallength=10;originalbedline=1;subregion=0,0	.	-
NC_010241.1	12	13	length=1;counts=3/1/2;originallength=10;originalbedline=1;subregion=2,2	.	-
NC_010241.1	14	15	length=1;counts=5/2/3;originallength=10;originalbedline=1;subregion=4,4	.	-
NC_010241.1	16	17	length=1;counts=10/8/2;originallength=10;originalbedline=1;subregion=6,6	.	-
NC_010241.1	18	19	length=1;counts=8/1/7;originallength=10;originalbedline=1;subregion=8,8	.	-
};

my $expected_output_min6 = qq{NC_010241.1	16	17	length=1;counts=10/8/2;originallength=10;originalbedline=1;subregion=6,6	.	-
NC_010241.1	18	19	length=1;counts=8/1/7;originallength=10;originalbedline=1;subregion=8,8	.	-
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

($return,$stdout,$stderr)=run_script('../scripts/049_bed2signal.pl' ,[$filename, 6]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With single file/condition is should exit with exit code 0');

$got	  = &bed::parser($stdout);
$expected = &bed::parser($expected_output_min6);
is_deeply($got,$expected,'output as expected');

done_testing();
