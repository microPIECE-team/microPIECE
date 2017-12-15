#!/usr/bin/env perl

use strict;
use warnings;

my %h = ();

while (<>)
{
    chomp;

    my @fields = split("\t", $_);

    next unless ($fields[2] =~ /mRNA|exon|CDS/);

    $fields[8] =~ /Dbxref=GeneID:(\d+).*Genbank:([^,;]+)/;

    my $name = $2;
    my $gene = $1;
    
    $fields[8] =~ /ID=([^;]+)/;
    my $id = $1;
    $fields[8] =~ /Parent=([^;]+)/;
    my $parent = $1;
    
    if ($fields[2] eq "mRNA")
    {
	my $mrna = $id;
	$h{$gene}{$mrna} = {};
    }
    elsif ($fields[2] eq "exon" && exists $h{$gene}{$parent} )
    {
	my $mrna = $parent;
	$h{$gene}{$mrna}{exon}{name} = $name;
	$h{$gene}{$mrna}{exon}{len} += abs($fields[3]-$fields[4]);
    }
    elsif (exists $h{$gene}{$parent})
    {
	my $mrna = $parent;
	$h{$gene}{$mrna}{cds}{name} = $name;
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
