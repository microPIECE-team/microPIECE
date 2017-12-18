#! /usr/bin/perl
use strict;
use warnings;

my $fasta_file	= $ARGV[0];
my $counter	= 0;
open(FA,"<",$fasta_file) || die;

my @header = ();
my $seq = "";

while(<FA>){
	chomp;
	my $fa_line	= $_;

	if($fa_line =~ /^>/ || eof(FA) )
	{
	    if (eof(FA) && $fa_line !~ /^>/)
	    {
		$seq .= $_;
	    }

	    # print one sequence per header
	    if ($seq)
	    {
		foreach my $head (@header)
		{
		    printf ">%s-%d %s\n%s\n", $head->{name}, ++$counter, $head->{location}, uc($seq);
		}
	    }

	    if ($fa_line =~ /annotation=([^-:;]+)/)
	    {
		my $annotationlist = $1;
		my $location = "unknown location"; # default location
		if ($fa_line =~ /::(.+:\d+-\d+\([-+]\))$/)
		{
		    $location = $1;
		}

		@header = map { {name => $_, location => $location} } (split(",", $annotationlist));

	    } elsif ($fa_line =~ /counts=([^-:;]+)/)
	    {
		print STDERR "No annotation found for line '$fa_line'! Skipping entry...\n";
		@header = (); # empty header array results in no output for that fasta block
	    } else {
		@header = ($fa_line);
	    }

	    $seq = "";

	} else {
	    $seq .= $_;
	}
}

close(FA) || die;
