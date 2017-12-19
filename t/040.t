#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

my($return,$stdout,$stderr)=run_script('../scripts/046_merge_bed_files.pl');
like($stderr,qr/Example defining to classes with two files each and output to merged.bed/,'print help when no args provided');

($return,$stdout,$stderr)=run_script('../scripts/046_merge_bed_files.pl',[qw(--input bed1=040/041_testset_01.bed)]);

is($return,0,'script with one bed file as input');
print "$return\n$stdout\n$stderr\n";

done_testing();
