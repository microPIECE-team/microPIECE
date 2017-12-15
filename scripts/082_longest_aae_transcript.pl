#!/usr/bin/env perl

use strict;
use warnings;

my %h = ();

while (<>)
{
    chomp;
    
    my @F = split("\t", $_);
    
    next unless ($F[2] =~ /mRNA|exon|CDS/);

    $F[8] =~ /Dbxref=VectorBase:([^-]+)-[RP]([^,]+).+Genbank:([^,]+)/;

    if ($F[2] eq "mRNA")
    { 
	$h{$1}{$2} = {}; 
    }
    elsif ($F[2] eq "exon" && exists $h{$1}{$2} ) 
    { 
	$h{$1}{$2}{exon}{name} = $3; 
	$h{$1}{$2}{exon}{len} += abs($F[3]-$F[4]); 
    } 
    elsif (exists $h{$1}{$2})
    { 
	$h{$1}{$2}{cds}{name} = $3; 
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
