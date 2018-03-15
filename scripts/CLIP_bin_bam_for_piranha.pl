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
	my $next_val = 0;
	for(my $start=1; $start<=$seq->{len}; $start+=$binsize)
	{
	    my $stop = $start+$binsize-1;
	    if ($stop>$seq->{len})
	    {
		$stop = $seq->{len};
	    }
	    if ($start>=$next_val)
	    {
		printf STDERR "\rWorking on %s:%d-%d (finished: %.1f%%)", $seq->{name}, $start, $stop, ($start/$seq->{len}*100);
		$next_val = $start+$seq->{len}/1000;
	    }
	    foreach my $strand ("-", "+")
	    {
		# the bin starts at $start and ends at $stop
		my $strand_samtools = "-F 0x10";
		if ($strand eq "+")
		{
		    $strand_samtools = "-f 0x10"
		}
		print join("\t", ($seq->{name}, $start-1, $stop, ".", ".", $strand)), "\n";
	    }
	}
	print STDERR " FINISHED\n";
    }
}

sub parse_header
{
    my ($header) = @_;

    my @seq = map { /SN:(\S+)\s+LN:(\d+)/; { name => $1, len => int($2) }; } (grep {/^\@SQ\s+SN:/} (split("\n", $$header)));
    @seq = sort { $a->{name} cmp $b->{name} } (@seq);

    return \@seq;
}

sub check4fastforward{
    my ($start, $fast_forward, $bam, $seq, $bin_size) = @_;

    while ($fast_forward>5*$binsize)
    {
	if ($start+$fast_forward>$seq->{len})
	{
	    return;
	}
	my $stop = $start+$fast_forward-1;

	my $cmd = sprintf("samtools view -c -F 0x4 %s %s:%d-%d", $bam, $seq->{name}, $start, $stop);
	my $counts =  qx($cmd)+0;
	if ($? != 0)
	{
	    die "Error running command '$cmd'\n";
	}

	return $fast_forward if ($counts == 0);
	$fast_forward = int($fast_forward>>1);
    }

    return;
}
