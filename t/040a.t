#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

use bed;

my $first_file = qq{NC_010241.1	0	1	read0	.	+	.
NC_010241.1	10	20	read1	.	+	.
NC_010241.1	10	20	read2	.	-	.
NC_010241.1	10	20	read3	.	+	.
NC_010241.1	30	40	read4	.	+	.
NC_010241.1	39	45	read5	.	+	.
NC_010241.1	150	160	read6	.	+	.
NC_010241.1	160	165	read7	.	+	.
NC_010241.1	70	80	read8	.	+	.
};

my $expected_output = qq{NC_010241.1	0	1	length=1;counts=1/1	.	+
NC_010241.1	30	45	length=15;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,2/2,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	10	20	length=10;counts=2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2	.	+
NC_010241.1	10	20	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	-
NC_010241.1	150	165	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	70	80	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+
};

# remove empty lines
$first_file =~ s/^\s*$//mg;
$expected_output =~ s/^\s*$//mg;

my ($fh, $filename) = File::Temp::tempfile();

print $fh $first_file;
close($fh) || die;

my($return,$stdout,$stderr)=run_script('../scripts/046_merge_bed_files.pl');
like($stderr,qr/Example defining to classes with two files each and output to merged.bed/,'print help when no args provided');

($return,$stdout,$stderr)=run_script('../scripts/046_merge_bed_files.pl',["--input", "bed1=$filename"]);

my $got	= &bed::parser($stdout);
my $expected = &bed::parser($expected_output);
is_deeply($got,$expected,'output as expected');

done_testing();
