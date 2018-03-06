#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

use bed;

my $first_file = qq{NC_010241.1	10	15	read1	.	+	.
NC_010241.1	50	55	read3	.	+	.
NC_010241.1	20	25	read5	.	+	.
};

my $second_file = qq{NC_010241.1	10	15	read2	.	-	.
NC_010241.1	45	51	read4	.	+	.
NC_010241.1	21	23	read6	.	+	.
};

my $expected_output = qq{NC_010241.1	10	15	length=5;counts=1/1/0,1/1/0,1/1/0,1/1/0,1/1/0	.	+
NC_010241.1	10	15	length=5;counts=1/0/1,1/0/1,1/0/1,1/0/1,1/0/1	.	-
NC_010241.1	45	55	length=10;counts=1/0/1,1/0/1,1/0/1,1/0/1,1/0/1,2/1/1,1/1/0,1/1/0,1/1/0,1/1/0	.	+
NC_010241.1	20	25	length=5;counts=1/1/0,2/1/1,2/1/1,1/1/0,1/1/0	.	+
};

my ($fh1, $filename1) = File::Temp::tempfile();

print $fh1 $first_file;
close($fh1) || die;

my ($fh2, $filename2) = File::Temp::tempfile();

print $fh2 $second_file;
close($fh2) || die;

my ($return,$stdout,$stderr)=run_script('../scripts/CLIP_merge_bed_files.pl',["--input", "bed1=$filename1", '--input', "bed2=$filename2"]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With multiple files for a single condition it should exit with exit code 0');

my $got	= &bed::parser($stdout);
my $expected = &bed::parser($expected_output);
is_deeply($got,$expected,'output as expected');

done_testing();
