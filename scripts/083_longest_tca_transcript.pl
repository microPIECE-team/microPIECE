#!/usr/bin/env perl

use strict;
use warnings;

my %h = ();

while (<>)
{
    chomp;

    my @fields = split("\t", $_);

    next unless ($fields[2] =~ /mRNA|exon|CDS/);

    $fields[8] =~ /ID=([^;]+).*Dbxref=GeneID:(\d+).*Genbank:([^,;]+)/;

    if ($fields[2] eq "mRNA")
    {
	$h{$2}{$1} = {};
    }
    elsif ($fields[2] eq "exon" && exists $h{$2}{$1} )
    {
	$h{$2}{$1}{exon}{name} = $3;
	$h{$2}{$1}{exon}{len} += abs($fields[3]-$fields[4]);
    }
    elsif (exists $h{$2}{$1})
    {
	$h{$2}{$1}{cds}{name} = $3;
    }
    else
    {
	warn "For line # $.('$_') no mRNA was found\n";
    }
}

foreach my $vec(sort keys %h)
{
    my ($longest) = sort {$h{$vec}{$b}{exon}{len} <=> $h{$vec}{$a}{exon}{len}} (keys %{$h{$vec}});

    next unless (exists $h{$vec}{$longest}{exon} && exists $h{$vec}{$longest}{cds});

    print join("\t", sort ($h{$vec}{$longest}{cds}{name}, $h{$vec}{$longest}{exon}{name})), "\n";
}
