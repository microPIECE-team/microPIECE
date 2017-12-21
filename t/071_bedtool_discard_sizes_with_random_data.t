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

my $input = generate_input(100);

my ($fh, $filename) = File::Temp::tempfile();

print $fh @{$input->{in}};
close($fh) || die;

my ($return,$stdout,$stderr)=run_script('../scripts/071_bedtool_discard_sizes.pl');
is(Test::Script::Run::last_script_exit_code(), 2, 'Without arguments is should exit with exit code 2');

($return,$stdout,$stderr)=run_script('../scripts/071_bedtool_discard_sizes.pl' ,[$filename, $min_val, $max_val]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With single file/condition is should exit with exit code 0');

my $got	= [ sort split("\n", $stdout) ];
my $expected = [ sort split("\n", join("", @{$input->{expected}})) ];
is_deeply($got,$expected,'output as expected');

done_testing();

sub generate_input
{
    my ($num) = @_;

    my %ret_val = ();

    for(my $i=0; $i<$num; $i++)
    {
	# two random integer numbers, first the length of the bed region (1-10000 bp),
	# 12.5 % below the minimum
	# 75 % inside the allowed range
	# 12.5 % above the maximum
	my $select = rand();
	my $max_len = 10000;
	my $len = 0;
	if ($select<0.125)
	{
	    $len = int(rand($min_val));
	} elsif ($select <=(0.75+0.125)) {
	    $len = int($min_val+rand($max_val-$min_val+1));
	} else {
	    $len = int($max_val+rand($max_len-$max_val+1));
	}

	# second a start (1-5000000)
	my $start = int(rand(5000000)+1);

	my $strand = '-';
	# 50 % probability to have + strand
	if (rand()>=0.5)
	{
	    $strand = "+";
	}

	my $bedline = join("\t", (".", $start, $len+$start-1, ".", ".", $strand))."\n";
	push(@{$ret_val{in}}, $bedline);

	if ($len >= $min_val && $len <= $max_val)
	{
	    push(@{$ret_val{expected}}, $bedline);
	}
    }

    return \%ret_val;
}
