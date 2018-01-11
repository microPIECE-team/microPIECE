#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

my %exp_fq	= (
"\@NB500930:5:H37JNAFXX:1:11101:8624:1029 1:N:0:ATCACG" => 0,
"CGCAAGAGAGTTCGTCGGGGCATGGAATTCTCGGGTGCCAAGGAACTCCAG" => 0,
"+" => 0,
"AAAAAEEEEEEEEEEEEEE/AEE6EEEAEEEEEAAEEEAEEEEEEE<<EEE" => 0,);




my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/071_filter_fastq_N.pl',[qw(070/071_test_in.fastq)]);

my %got = %{&parser($stdout)};
is_deeply (\%got,\%exp_fq,'071 output as expected');

#my $expected = &parser(join('',<DATA>));
#is_deeply($got,$expected,'071 output as expected');

done_testing();





sub parser{
	my $p_stdout		=	$_[0];
	my $p_header		= 	"";
	my $p_counter		= 	0;
	my %p_hash;
	foreach my $fq_line ( split("\n",$p_stdout)){
		$p_hash{$fq_line} = 0;
	}
	return \%p_hash;
}

__DATA__
@NB500930:5:H37JNAFXX:1:11101:8624:1029 1:N:0:ATCACG
CGCAAGAGAGTTCGTCGGGGCATGGAATTCTCGGGTGCCAAGGAACTCCAG
+
AAAAAEEEEEEEEEEEEEE/AEE6EEEAEEEEEAAEEEAEEEEEEE<<EEE
