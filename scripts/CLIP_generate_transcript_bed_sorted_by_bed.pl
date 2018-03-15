#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $gff             = undef;
my $sortby          = undef;
my $feature         = "exon";
my $progressenabled = 1;

GetOptions(
    'g|gff=s'     => \$gff,
    's|sortby=s'  => \$sortby,
    'f|feature=s' => \$feature,
    'p|progress!' => \$progressenabled,
    ) || die;

if ($progressenabled)
{
    require Term::ProgressBar;
}

die "Missing option for gff input file (specify via --gff FILENAME)\n" unless defined $gff;
die "Missing option for sort input file (specify via --sortby FILENAME)\n" unless defined $sortby;

die "Filename provided via --gff option is not accessable (provided name is '$gff')\n" unless -e $gff;
die "Filename provided via --sortby option is not accessable (provided name is '$sortby')\n" unless -e $sortby;

# parse the sort order
my %seen  = ();
my @order = ();
my $progress;
if ($progressenabled)
{
    $progress = Term::ProgressBar->new({name => "Order import", count => -s $sortby, remove => 0, ETA => 'linear'});
    $progress->minor(0);
}
my $next_update = 0;
open(FH, "<", $sortby) || die "Unable to open file '$sortby': $!\n";
while(<FH>)
{
    my $pos = tell(FH);
    if ($progressenabled)
    {
	$next_update = $progress->update($pos) if ($pos >= $next_update);
    }

    next if (/^#/);
    # a line has to be splittable at space/tab
    next unless (/^(\S+)\s*/);

    my $chr = $1;

    if (! exists $seen{$chr})
    {
	push(@order, $chr);
	$seen{$chr}{pos} = int(@order);
    }
    $seen{$chr}{counter}++;
}
close(FH) || die "Unable to close file '$sortby': $!\n";
if ($progressenabled)
{
    $progress->update(-s $sortby) if ((-s $sortby)>=$next_update);
}

# print a status
print STDERR "Position\tName\t# seen\n";
foreach my $item (sort { $seen{$a}{pos} <=> $seen{$b}{pos} } (keys %seen))
{
    printf STDERR "%d\t%s\t%d\n", $seen{$item}{pos}, $item, $seen{$item}{counter};
}

# open the gff and run through the lines, skipping everything, but the
# requested feature (exon by default)
if ($progressenabled)
{
    $progress = Term::ProgressBar->new({name => "Import of GFF", count => -s $gff, remove => 0, ETA => 'linear'});
    $progress->minor(0);
}
$next_update = 0;
my @feature4bed = ();
open(FH, "<", $gff) || die "Unable to open gff input file: $!\n";
while(<FH>)
{
    my $pos = tell(FH);
    if ($progressenabled)
    {
	$next_update = $progress->update($pos) if ($pos >= $next_update);
    }

    next if (/^#/);
    my ($chr, undef, $feat, $start, $stop, undef, $strand) = split(/\t/, $_);

    next unless ($feat eq $feature);

    die "Chromosome '$chr' was not in sortby file\n" unless (exists $seen{$chr});

    push(@feature4bed, {
	order => $seen{$chr}{pos},
	str   => join("\t", ($chr, $start-1, $stop, ".", ".", $strand)),
	start => $start,
	 } );
}
close(FH) || die "Unable to close gff input file: $!\n";

if ($progressenabled)
{
    $progress->update(-s $gff) if ((-s $gff)>=$next_update);
}

# Sorting the input according to the sort by file
foreach my $entry (sort { $a->{order} <=> $b->{order} || $a->{start} <=> $b->{start} || $a->{str} cmp $b->{str} } (@feature4bed))
{
    print $entry->{str}, "\n";
}
