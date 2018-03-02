#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

my %cluster = ();

while(<>)
{
    next if (/^#/);
    
    chomp;
    
    my @fields = split("\t", $_);

    my $data = {};

    next unless (defined $fields[3]);
    
    foreach my $current (split(";", $fields[3]))
    {
	my ($name, $value) = split("=", $current, 2);

	$data->{$name} = $value;
    }

    # split the counts
    die "Missing counts" unless (exists $data->{counts});

    $data->{counts} = [ split(",", $data->{counts}) ];

    my $uniq_pattern = {};
    foreach my $pattern (@{$data->{counts}})
    {
	$uniq_pattern->{$pattern}++;
    }

    my $uniq_cluster = {};
    
    foreach my $pattern (keys %{$uniq_pattern})
    {
	my @conditions = split("/", $pattern);

	for(my $i=1; $i<@conditions; $i++)
	{
	    next unless ($conditions[$i] > 0);
	    $uniq_cluster->{$i}++;
	}
    }
    foreach my $cluster (keys %{$uniq_cluster})
    {
	push(@{$cluster{$cluster}}, sprintf("line%05d", $.));
    }
}

foreach my $condition (sort keys %cluster)
{
    printf "%s:%s;\n", sprintf("condition%d", $condition), join(",", @{$cluster{$condition}});
}
