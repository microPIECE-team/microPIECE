#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

my $expected_data	= qq{\@NB500930:5:H37JNAFXX:1:11101:8624:1029 1:N:0:ATCACG
CGCAAGAGAGTTCGTCGGGGCATGGAATTCTCGGGTGCCAAGGAACTCCAG
+
AAAAAEEEEEEEEEEEEEE/AEE6EEEAEEEEEAAEEEAEEEEEEE<<EEE};




my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/071_filter_fastq_N.pl',["070/071_test_in.fastq"]);
isnt($stdout,"","071 output is not empty");

my $got = join("\n",sort(split("\n",$stdout)));
my $exp = join("\n",sort(split("\n",$expected_data)));

is($got,$exp,'071 output as expected');


done_testing();





