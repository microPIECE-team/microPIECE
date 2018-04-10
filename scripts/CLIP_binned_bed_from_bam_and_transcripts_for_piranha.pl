#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;

=pod

=head1 NAME

CLIP_binned_bed_from_bam_and_transcripts_for_piranha.pl - Pre binning of a BAM file with optional pseudocounts for transcript sites

=head1 DESCRIPTION

This script is used to divide a BAM file into bins of a specified size
(default 20). Additionally, pseudocounts can be added on every
position covered by a transcript. Therefore, an input GFF file can be
specified. Output is printed to STDOUT.

=head1 SYNOPSIS

   CLIP_binned_bed_from_bam_and_transcripts_for_piranha.pl [options] --bam bamfile

   Options:
     --help
     --bam
     --size
     --transcripts
     --reqfeature

   # Example defining to classes with two files each and output to merged.bed
   CLIP_binned_bed_from_bam_and_transcripts_for_piranha.pl
      --bam         INPUT-BAM.bam \
      --size        32 \
      --transcripts ANNOTATION.gff \
      --reqfeature  mRNA

=head1 OPTIONS

=over 8

=item C<--bam|-b>

Specifies the input file. Should be a BAM file with a prepared index. Use
C<samtools index BAMFILENAME> to generate the required index.
Specifying a bam is mandatory.

=item C<--size|-s>

Specifies the binsize. A bins counter is increased for each read start
located inside that bin. This is following the binning/counting of
C<Piranha>. Bins with a count size of C<0> are excluded from the
output.

=item C<--transcripts|-t>

Specifies the transcript annotations. It has to be a GFF3
file. Required fields are C<chr (column 1)>, C<feature (column 3)>,
C<start (column 4)>, C<stop (column 5)>, and C<strand (column
7)>. Other fields are ignored.  Additionally, only lines with a
feature matching the C<--reqfeature> (default exon) are considered.  A
pseudocount of 1 is added to all bins, containing one of required
feature.

=item C<--reqfeature|-r>

Specifies the feature required from the transcript annotations. The
default value is "exon", meaning only lines with exon as feature will
result in a pseudocount of 1 for all bins covered.

=item C<--help|-h|-?>

Print that help message

=item C<--version|-V>

Print that help message

=back

=cut

use Getopt::Long;

use version 0.77; our $VERSION = version->declare("v1.5.0");
my $version      = 0;
my $help         = 0;

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
    'V|version'        => \$version,
    'h|help|?'         => \$help,
    ) || pod2usage(2);

pod2usage(1) if $help;
if ($version)
{
    print $VERSION, "\n";
    exit;
}

# definde the bams and transcripts if comma seperated
@bamfiles = split(",", join(",", @bamfiles));
if (@bamfiles>1)
{
    warn("Multiple bam files as input is currently not supported. Will create bin only for first bam: ".$bamfiles[0]." and skip the following: ".join(",", @bamfiles[1..@bamfiles-1])."\n");
    @bamfiles = ($bamfiles[0]);
}
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
    my $progress = Term::ProgressBar->new({name => "Transcript import from '$transcript'", count => -s $transcript, remove => 0, ETA => 'linear'});
    $progress->minor(0);
    my $next_val = 0;

    open(FH, "<", $transcript) || die "Unable to open transcript file '$transcript'\n";
    while (<FH>)
    {
	my $pos = tell(FH);
	$next_val = $progress->update($pos) if ($pos >= $next_val);
	chomp($_);

	# skip header
	next if (/^#/);

	# as the input is a gff file, we want to consider only the required feature
	my ($chr, undef, $feature, $start, $end, undef, $strand) = split(/\t/, $_);
	next unless (defined $feature && $feature eq $reqfeature);

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
	push(@{$imported_transcripts{$chr}{$strand}}, ($start, $end));

	$pseudocounts = 1;
    }
    close(FH) || die "Unable to close transcript file '$transcript'\n";
    #$progress->update(-s $transcript) if (-s $transcript <= $next_val);
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
	my $progress = Term::ProgressBar->new({name => $seq->{name}, count => $seq->{len}, remove => 0, ETA => 'linear'});
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

	    # set the progressbar
	    $next_val = $progress->update($start) if ($start >= $next_val);

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

	    $counts[$bin][$strand]++;
	}
	close(FH) || die "Unable to close pipe to '$cmd': $!\n";
	$progress->update($seq->{len}) if ($seq->{len} <= $next_val);

	if ($pseudocounts)
	{
	    my $list_of_feature_bins = expand_transcripts(\%imported_transcripts, $seq->{name});
	    my $bins2change = 0;
	    foreach my $strand (keys %{$imported_transcripts{$seq->{name}}})
	    {
		$bins2change += @{$list_of_feature_bins->[$strand]}
	    }
	    my $progress = Term::ProgressBar->new({name => $seq->{name}." Pseudocounts", count => $bins2change, remove => 0, ETA => 'linear'});
	    $progress->minor(0);
	    my $next_val = 0;
	    my $bins_changed = 0;
	    foreach my $strand (keys %{$imported_transcripts{$seq->{name}}})
	    {
		foreach my $bin (@{$list_of_feature_bins->[$strand]})
		{
		    $counts[$bin][$strand]++;
		    $bins_changed++;
		    $next_val = $progress->update($bins_changed) if ($bins_changed >= $next_val);
		}
	    }
	   # $progress->update($bins2change) if ($bins2change <= $next_val);
	}

	$progress = Term::ProgressBar->new({name => $seq->{name}." BEDoutput", count => int(@counts), remove => 0, ETA => 'linear'});
	$progress->minor(0);
	$next_val = 0;
	for(my $i=0; $i<@counts; $i++)
	{
	    $next_val = $progress->update($i) if ($i >= $next_val);

	    # get the start and end coordinate
	    my $start = $i*$binsize;
	    my $stop = $start+$binsize;

	    foreach my $strand (0, 1)
	    {
		if (defined $counts[$i] && defined $counts[$i][$strand] && $counts[$i][$strand] > 0)
		{
		    print join("\t", ($seq->{name}, $start, $stop, ".", $counts[$i][$strand], ($strand == 0) ? "+" : "-")), "\n";
		}
	    }
	}
	$progress->update(@counts+0) if (@counts <= $next_val);
    }
}

sub parse_header
{
    my ($header) = @_;

    my @seq = map { /SN:(\S+)\s+LN:(\d+)/; { name => $1, len => int($2) }; } (grep {/^\@SQ\s+SN:/} (split("\n", $$header)));

    return \@seq;
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

sub expand_transcripts
{
    my ($transcripts, $seqname) = @_;
    my %seen = ();

    my @strands = keys %{$transcripts->{$seqname}};

    foreach my $strand (@strands)
    {
	for(my $i=0; $i<@{$transcripts->{$seqname}{$strand}}; $i+=2)
	{
	    my $start = $transcripts->{$seqname}{$strand}[$i];
	    my $stop  = $transcripts->{$seqname}{$strand}[$i+1];
	    for(my $bin=$start; $bin<=$stop; $bin++)
	    {
		$seen{$strand}{$bin}++;
	    }
	}
    }

    my @result = ();
    foreach my $strand (@strands)
    {
	$result[$strand] = [ keys %{$seen{$strand}} ];
    }

    return \@result;

}
