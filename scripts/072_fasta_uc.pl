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
		    printf ">%s-%d\n%s\n", $head, ++$counter, uc($seq);
		}
	    }

	    if ($fa_line =~ /annotation=([^-:;]+)/)
	    {
		@header = split(",", $1);
	    } elsif ($fa_line =~ /counts=([^-:;]+)/)
	    {
		print STDERR "No annotation found for line '$fa_line'! Skipping entry...\n";
		@header = ();
	    } else {
		@header = ($fa_line);
	    }

	    $seq = "";

	} else {
	    $seq .= $_;
	}
}

close(FA) || die;
