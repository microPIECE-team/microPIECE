#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $gff            = undef;
my $sortby         = undef;
my $feature        = "exon";

GetOptions(
    'g|gff=s'     => \$gff,
    's|sortby=s'  => \$sortby,
    'f|feature=s' => \$feature,
    ) || die;

die "Missing option for gff input file (specify via --gff FILENAME)\n" unless defined $gff;
die "Missing option for sort_order_bed input file (specify via --sortby FILENAME)\n" unless defined $sortby;

die "Filename provided via --gff option is not accessable (provided name is '$gff')\n" unless -e $gff;
die "Filename provided via --sortby option is not accessable (provided name is '$sortby')\n" unless -e $sortby;

# parse the sort order
my %seen  = ();
my @order = ();
open(FH, "<", $sort_order_bed) || die "Unable to open file '$sort_order_bed': $!\n";
while(<FH>)
{
    next if (/^#/);
    # a line has to be splittable at space/tab
    next unless (/^(\S+)\s/);

    if (! exists $seen{$next})
    {
	push(@order, $next);
	$seen{$next}{pos} = int(@order);
    }
    $seen{$next}{counter}++;
}
close(FH) || die "Unable to close file '$sort_order_bed': $!\n";

# print a status
print STDERR "Position\tName\t# seen\n";
foreach my $item (sort { $seen{$a}{pos} <=> $seen{$b}{pos} } (keys %seen))
{
    printf STDERR "%d\t%s\t%d\n", $seen{$item}{pos}, $item, $seen{$item}{counter};
}
