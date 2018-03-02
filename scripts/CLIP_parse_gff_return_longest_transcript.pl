#!/usr/bin/env perl

use strict;
use warnings;

# due to we use gffread to convert gff to proteins/transcripts, I
# assume, that each mRNA has an unique ID and that ID is used as
# parent for the CDS/exon

my %dat = ();
my %mrna2gene = ();

while (<>)
{
    chomp;

    next if ($_ =~ /^#/);

    my @fields = split("\t", $_);

    next unless ($fields[2] =~ /mRNA|exon/);

    if ($fields[2] =~ /mRNA/)
    {
	die unless ($fields[8] =~ /ID=([^;]+)/);
	my $id = $1;

	die unless ($fields[8] =~ /Parent=([^;]+)/);
	my $parent = $1;

	$dat{$parent}{$id} = 0;
	$mrna2gene{$id} = $parent;
    } else {
	die unless ($fields[8] =~ /Parent=([^;]+)/);
	my $mrna = $1;
	$dat{$mrna2gene{$mrna}}{$mrna} += abs($fields[3]-$fields[4])
    }
}

printf STDERR "Found %d mRNAs\n", int(keys %mrna2gene);

foreach my $gene (keys %dat)
{
    my @mrnas = sort {$dat{$gene}{$b} <=> $dat{$gene}{$b}} (keys %{$dat{$gene}});
    my $longest = shift @mrnas;

    print join("\t", ($longest, $longest)), "\n";
}
