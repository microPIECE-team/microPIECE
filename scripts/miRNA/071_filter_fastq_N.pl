#! /usr/bin/perl
use strict;
use warnings;
use GetOpt::Long;

# input is a FASTQ file, which will be scanned for sequences containing Ns.
# All FASTQ blocks containing Ns will be skipped. All other are printed to STDOUT.

while (<>)
{
	my $header  = $_;
	my $seq     = <> || die "wrong number of lines\n";
	my $header2 = <> || die "wrong number of lines\n";
	my $qual    = <> || die "wrong number of lines\n";

	# check for N or n in a string
	my $num_Ns = $seq =~ tr/Nn/Nn/;
	
	next if ($num_Ns > 0);
	
	print $header, $seq, $header2, $qual;
}
