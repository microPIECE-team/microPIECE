#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my @bamfiles = ();
my $binsize  = 20;

GetOptions(
    'b|bam=s@' => \@bamfiles,
    's|size=i' => \$binsize,
    ) || die;

# definde the bams if comma seperated
@bamfiles = split(",", join(",", @bamfiles));

# try to call samtools
my $samtools = qx(samtools --version);
if ($? != 0)
{
    die "Wrong version or no samtools available\n";
}

foreach my $bam (@bamfiles)
{
    my $cmd = "samtools view -H $bam";
    my $header = qx($cmd);
    if ($? != 0)
    {
	die "Error running command '$cmd'\n";
    }
    my $seq = parse_header(\$header);

    foreach my $seq (@{$seq})
    {
	print STDERR "Working on $seq->{name}:";
	# generate a bin of size binsize
	for(my $start=1; $start<=$seq->{len}; $start+=$binsize)
	{
	    my $stop = $start+$binsize-1;
	    if ($stop>$seq->{len})
	    {
		$stop = $seq->{len};
	    }
	    print STDERR "\rWorking on $seq->{name}:$start-$stop";
	    foreach my $strand ("-", "+")
	    {
		# the bin starts at $start and ends at $stop
		my $strand_samtools = "-F 0x10";
		if ($strand eq "+")
		{
		    $strand_samtools = "-f 0x10"
		}
		my $cmd = sprintf("samtools view -c -F 0x4 %s %s %s:%d-%d", $strand_samtools, $bam, $seq->{name}, $start, $stop);
		my $counts =  qx($cmd)+0;
		if ($? != 0)
		{
		    die "Error running command '$cmd'\n";
		}
		if ($counts > 0)
		{
		    # print BED line

		    # Start needs to be decreased, due to 0-based BED
		    # vs. 1-based BAM coordinates, but $stop requires no
		    # modification, due to it is exclusive in BED
		    print join("\t", ($seq->{name}, $start-1, $stop, ".", $counts, $strand)), "\n";
		}
	    }
	}
	print STDERR "FINISH\n";
    }
}

sub parse_header
{
    my ($header) = @_;

    my @seq = map { /SN:(\S+)\s+LN:(\d+)/; { name => $1, len => int($2) }; } (grep {/^\@SQ\s+SN:/} (split("\n", $$header)));
    @seq = sort { $a->{name} cmp $b->{name} } (@seq);

    return \@seq;
}
