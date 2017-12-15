#! /usr/bin/perl
use strict;
use warnings;

my $bedfile	= $ARGV[0];
my $min		= $ARGV[1];	# eg 22
my $max		= $ARGV[2];	# eg 50

open(BED ,"<",$bedfile) || die;
while(<BED>){
	chomp;

	# ignore comment line
	if ($_ =~ /^#/)
	{
	    print $_, "\n";
	    next;
	}

	my $bed_line	= $_;
	my @bed_array	= split(" ",$bed_line);
	my $bed_start	= $bed_array[1];
	my $bed_stop	= $bed_array[2];
	my $bed_len	= $bed_stop-$bed_start;
	next if($bed_len < $min);
	next if($bed_len > $max);
	print "$bed_line\n";
}
close(BED)|| die;


