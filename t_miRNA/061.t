#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
use Switch;
# result file ends with an empty line

# generate datastructure for expected data 
my %exp_rpm = (
"rpm;condition;miRNA"		=> 0,
"200000;a;tca-miR-3811c-1-5p"	=> 0,
"225000;a;tca-miR-1-3p"		=> 0,
"275000;a;tca-miR-1-5p"		=> 0,
"300000;a;tca-miR-3811c-2-5p"	=> 0,
"200000;b;tca-miR-3811c-1-5p"	=> 0,
"200000;b;tca-miR-1-3p"		=> 0,
"200000;b;tca-miR-1-5p"		=> 0,
"400000;b;tca-miR-3811c-2-5p"	=> 0);



my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/061_sam2de.pl',[qw(060/061_test.cfg 060/061_mature.fa)]);

my %got_rpm	= %{&parser($stdout)};

is_deeply(\%got_rpm,\%exp_rpm,'061 output as expected');
done_testing();



sub parser{
	my $p_stdout	= $_[0];
	my %p_rpm;
	foreach my $p_line (split("\n",$p_stdout)){
		$p_rpm{$p_line} = 0;
	}
	return(\%p_rpm);
}


