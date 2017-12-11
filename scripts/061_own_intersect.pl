#!/usr/bin/env perl

use strict;
use warnings;

my %input = ();

use Getopt::Long;
use Data::Dumper;

GetOptions(
    "input=s%" => \%input
    );

# split input keys into seperate files
foreach my $key (keys %input)
{
    $input{$key} = [ split(/,/, $input{$key}) ];
}

print Dumper(\%input);
