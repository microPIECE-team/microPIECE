#!/usr/bin/env perl

use strict;
use warnings;

my %parent2childs = ();
my %id2index = ();
my @dat = ();

while (<>)
{
    chomp;

    next if ($_ =~ /^#/);

    my @fields = split("\t", $_);

    next unless ($fields[2] =~ /mRNA|exon|CDS|gene/);

    push(@dat, {orig => \@fields, type => $fields[2]});

    my $index = @dat-1;

    # get id
    $dat[$index]{orig}[8] =~ /ID=([^;]+)/;
    my $id = $1;
    $dat[$index]{id} = $id;
    # get parent
    $dat[$index]{orig}[8] =~ /Parent=([^;]+)/;
    my $parent = $1;
    $dat[$index]{parent} = $parent;

    $id2index{$id} = $index;
    push(@{$parent2childs{$parent}}, $id);

    unless (exists $id2index{$parent})
    {
	warn "Missing parent '$parent': line: $.\n";
	next;
    }

    if ($fields[2] eq "exon")
    {

	$dat[$id2index{$parent}]{len} += abs($fields[3]-$fields[4]);
	push(@{$dat[$id2index{$parent}]{exons}}, $dat[$index]);
    }elsif ($fields[2] eq "CDS")
    {
	push(@{$dat[$id2index{$parent}]{cds}}, $dat[$index]);
    }
}

my @mrnas = map { $_->{id} } grep {$_->{type} eq "mRNA"} (@dat);

printf STDERR "Found %d mRNAs\n", @mrnas+0;

my %genes = ();

for(my $i=0; $i<@mrnas; $i++)
{
    my $index4mrna = $id2index{$mrnas[$i]};
    my $parent4mrna = $dat[$index4mrna]{parent};

    push(@{$genes{$parent4mrna}}, $index4mrna);
}

printf STDERR "Found %d genes\n", int(keys %genes);

foreach my $gene (sort keys %genes)
{
    # printf STDERR "Gene: %s\n", $gene;
    my $ref_mrnas = [ grep { $dat[$id2index{$_}]{type} eq "mRNA" } (@{$parent2childs{$gene}}) ];

    my @mrnas_sorted_by_length = map { $dat[$_] } sort { $dat[$b]{len} <=> $dat[$a]{len} || $dat[$b]{id} cmp $dat[$a]{id} } map {$id2index{$_}} (@{$ref_mrnas});

    # build structure:
    my $struct = { gene => $gene, mrnas => [] };

    foreach my $mrna (@mrnas_sorted_by_length)
    {
	my $cds = [ map {$_->{id} } (@{$mrna->{cds}}) ];
	my $exons = [ map {$_->{id} } (@{$mrna->{exons}}) ];

	push(@{$struct->{mrnas}}, { cds => $cds, exons => $exons });
    }

    my $longest_mrna = $mrnas_sorted_by_length[0];

    my ($protid, $mrnaid);

    # get Protein-ID from cds
    if ($longest_mrna->{cds}[0]{orig}[8] && $longest_mrna->{cds}[0]{orig}[8] =~ /Genbank:([^,;]+)/)
    {
	$protid = $1;
    } else {
	print STDERR "Something wrong with proteinid in ".Dumper($longest_mrna->{cds}[0]);
    }
    # get mRNA-ID from exons
    if ($longest_mrna->{exons}[0]{orig}[8] && $longest_mrna->{exons}[0]{orig}[8] =~ /Genbank:([^,;]+)/)
    {
	$mrnaid = $1;
    } else {
	print STDERR "Something wrong with mrnaid in ".Dumper($longest_mrna->{exons}[0]);
    }

    if ($mrnaid && $protid)
    {
	print join("\t", $mrnaid, $protid), "\n";
    }
}

use Data::Dumper;

# foreach my $vec(sort keys %h)
# {
#     my ($longest) = sort {$h{$vec}{$b}{exon}{len} <=> $h{$vec}{$a}{exon}{len}} (keys %{$h{$vec}});

#     next unless (exists $h{$vec}{$longest}{exon} && exists $h{$vec}{$longest}{cds});

#     print join("\t", sort ($h{$vec}{$longest}{cds}{name}, $h{$vec}{$longest}{exon}{name})), "\n";
# }
