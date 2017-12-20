#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

# the script only works on fields 2&3 of the bed file, therefore all
# other information are irrelevant and should be random values

# set a fixed value for the prng
srand(29164);

my $number_of_lines = 100;
my $min_val = 22;
my $max_val = 50;

my $expected_output = ();

my $input = generate_input(100);

print Dumper($input); use Data::Dumper;

sub generate_input
{
    my ($num) = @_;

    my @ret_val = ();

    for(my $i=0; $i<$num; $i++)
    {
	# two random integer numbers, first the length of the bed region (1-10000 bp)
	my $len = int(rand(10000)+1);
	# second a start (1-5000000)
	my $start = int(rand(5000000)+1);

	my $strand = '-';
	# 50 % probability to have + strand
	if (rand()>=0.5)
	{
	    $strand = "+";
	}

	push(@ret_val, join("\t", (".", $start, $len+$start-1, ".", ".", $strand))."\n");
    }

    return \@ret_val;
}

my ($fh, $filename) = File::Temp::tempfile();

print @{$input};

print $fh @{$input};
close($fh) || die;

my ($return,$stdout,$stderr)=run_script('../scripts/071_bedtool_discard_sizes.pl');
is(Test::Script::Run::last_script_exit_code(), 2, 'Without arguments is should exit with exit code 2');

($return,$stdout,$stderr)=run_script('../scripts/071_bedtool_discard_sizes.pl' ,[$filename, $min_val, $max_val]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With single file/condition is should exit with exit code 0');

my $got	= [ split("\n", $stdout) ];
my $expected = [ split("\n", $expected_output) ];
is_deeply($got,$expected,'output as expected');

done_testing();
