#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my @bamfiles     = ();
my @transcripts  = ();
my $binsize      = 20;
my $reqfeature   = "exon";
my $pseudocounts = 0;  

GetOptions(
    'b|bam=s@'         => \@bamfiles,
    's|size=i'         => \$binsize,
    't|transcripts=s@' => \@transcripts,
    'r|reqfeature=s'   => \$reqfeature,
    ) || die;

# definde the bams and transcripts if comma seperated
@bamfiles = split(",", join(",", @bamfiles));
@transcripts = split(",", join(",", @transcripts));

# try to call samtools
my $samtools = qx(samtools --version);
if ($? != 0)
{
    die "Wrong version or no samtools available\n";
}

use Term::ProgressBar;

# import the transcripts
my %imported_transcripts = ();
foreach my $transcript (@transcripts)
{
    open(FH, "<", $transcript) || die "Unable to open transcript file '$transcript'\n";
    while (<FH>)
    {
	chomp($_);

	# as the input is a gff file, we want to consider only the required feature
	my ($chr, undef, $feature, $start, $end, undef, $strand) = split(/\t/, $_);
	next unless ($feature eq $reqfeature);

	# transform the start and end coordinates into BED format (0-based, end-exclusive)
	($start, $end) = transform2bed($start, $end);

	# strand + => 0; - => 1;
	if ($strand eq "+")
	{
	    $strand = 0;
	} elsif ($strand eq "-")
	{
	    $strand = 1;
	}

	# start and stop are transformed into bins
	$start = int($start/$binsize);
	$end   = int(($end-1)/$binsize);
	push(@{$imported_transcripts{$chr}{$strand}}, { start => $start, stop => $end });

	$pseudocounts = 1;
    }
    close(FH) || die "Unable to close transcript file '$transcript'\n";
}
# finally sort the imported transcripts by coordinate
foreach my $chr (keys %imported_transcripts)
{
    foreach my $strand (keys %{$imported_transcripts{$chr}})
    {
	my @bins_with_feature = ();
	
	foreach my $exon (@{$imported_transcripts{$chr}{$strand}})
	{
	    foreach my $bin ($exon->{start}..$exon->{stop})
	    {
		$bins_with_feature[$bin]++;
	    }
	}
	$imported_transcripts{$chr}{$strand} = \@bins_with_feature;
    }
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
    my $total_length = 0;
    foreach my $seq (@{$seq})
    {
	$total_length += $seq->{len};
    }

    foreach my $seq (@{$seq})
    {
	my $progress = Term::ProgressBar->new({name => $seq->{name}, count => $seq->{len}, remove => 0});
	$progress->minor(0);
	my $next_val = 0;

	my @counts = ();
	$#counts = int($seq->{len}/$binsize)+1;

	my $cmd = sprintf("samtools view -F 0x4 %s %s:%d-%d", $bam, $seq->{name}, 1, $seq->{len});
	open(FH, "$cmd |") || die "Unable to open pipe to '$cmd': $!\n";
	while (<FH>)
	{
	    # ignore lines starting with @
	    next if (/^@/);
	    
	    # got a mapped read
	    my ($name, $flag, $chr, $start, $qual) = split("\t", $_);

	    # is the read a forward mapped read
	    # assume forward strand
	    my $strand = 0;
	    if ($flag & 16)
	    {
		# mapping is reverse
		$strand = 1;
	    }

	    ($start, undef) = transform2bed($start);

	    # according to Piranha the start position determines the bin
	    my $bin = int($start/$binsize);
	    
	    # set the progressbar
	    $next_val = $progress->update($start) if ($start >= $next_val);

	    $counts[$bin][$strand]++;
	}
	close(FH) || die "Unable to close pipe to '$cmd': $!\n";

	if ($pseudocounts)
	{
	    foreach my $strand (keys %{$imported_transcripts{$seq->{name}}})
	    {
		foreach my $bin (@{$imported_transcripts{$seq->{name}}{$strand}})
		{
		    $counts[$bin][$strand]++;
		}
	    }
	}
	
	for(my $i=0; $i<@counts; $i++)
	{
	    # get the start and end coordinate
	    my $start = $i*$binsize;
	    my $stop = $start+$binsize;

	    if (defined $counts[$i][0] && $counts[$i][0] > 0)
	    {
		print join("\t", ($seq->{name}, $start, $stop, ".", $counts[$i][0], "+")), "\n";
	    }

	    if (defined $counts[$i][1] && $counts[$i][1] > 0)
	    {
		print join("\t", ($seq->{name}, $start, $stop, ".", $counts[$i][0], "-")), "\n";
	    }

	}
    }
}

sub parse_header
{
    my ($header) = @_;

    my @seq = map { /SN:(\S+)\s+LN:(\d+)/; { name => $1, len => int($2) }; } (grep {/^\@SQ\s+SN:/} (split("\n", $$header)));

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

sub transform2bed
{
    my ($start, $end) = @_;

    if (defined $start)
    {
	$start--;
    }
    
    if (defined $end)
    {
	$end = $end;
    }

    return($start, $end);
}

sub transformfrombed
{
    my ($start, $end) = @_;

    if (defined $start)
    {
	$start++;
    }
    
    if (defined $end)
    {
	$end = $end;
    }

    return($start, $end);
}
