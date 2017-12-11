#!/usr/bin/env perl

use strict;
use warnings;

my %input = ();

use Getopt::Long;
use Data::Dumper;

GetOptions(
    "input=s%" => \%input
    );

# split input keys into seperate files
foreach my $key (keys %input)
{
    $input{$key} = [ split(/,/, $input{$key}) ];
}

my %genome = ();

foreach my $key (keys %input)
{
    foreach my $file (@{$input{$key}})
    {
	warn "Working on file '$file'\n";
	open(FH, "<", $file) || die "Unable to open file '$file': $!";
	while(<FH>)
	{
	    # go through the bed files
	    chomp;

	    # NW_001809801.1	79510	79534	XM_001647792.1	1	+
	    # NW_001809801.1	79527	79534	XM_001647792.1	1	+
	    # NW_001809801.1	248783	248807	XM_001647796.1	1	-
            # NW_001809801.1	533942	533950	XM_001647802.1	1	-
            # NW_001809801.1	2779750	2779790	XM_001647839.1	1	+

	    my ($chromosome, $start, $stop, $mrna_ids, undef, $strand) = split(/\t/, $_);

	    # for each chromosomal position we need to consider:
	    # Counter + strand
	    # Counter - strand

	    for (my $i=$start; $i<=$stop; $i++)
	    {
		$genome{$chromosome}[$i]{$strand}++;
	    }
	}
	close(FH) || die "Unable to close file '$file': $!";
    }
}

# print the output
foreach my $chromosome (sort keys %genome)
{
    foreach my $strand (qw(- +))
    {
	my $start = -1;
	for (my $i=0; $i<@{$genome{$chromosome}}; $i++)
	{
	    if (defined $genome{$chromosome}[$i] && exists $genome{$chromosome}[$i]{$strand} && $genome{$chromosome}[$i]{$strand}>0)
	    {
		# we found a new block
		$start = $i;
		my $stop = -1;
		for (my $j=$i+1; $j<@{$genome{$chromosome}}; $j++)
		{
		    if (defined $genome{$chromosome}[$j] && exists $genome{$chromosome}[$j]{$strand} && $genome{$chromosome}[$j]{$strand}>0)
		    {
			$stop = $j;
		    } else {
			last;
		    }
		}
		if ($stop == -1)
		{
		    $stop = int(@{$genome{$chromosome}});
		}

		print join("\t", ($chromosome, $start, $stop, ".", ".", $strand)), "\n";
		$i = $stop+1;
		$start = -1; $stop = -1;
	    }
	}
    }
}

#print Dumper(\%genome);
