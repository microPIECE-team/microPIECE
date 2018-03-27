package mining;

use strict;
use warnings;

use Digest::MD5;

=pod

=head2 SUBROUTINE C<parse_mirbase_dat()>

Will parse the *.dat files from mirBASE and return a reduced structure.

=head3 PARAMETER

=over 4

=item C<infile>

A filename of the input file

=item C<species>

A three letter code for a species. If one is provided, only species
matching the three letter code are returned. This filter step is
skipped, if species is not provided (or C<undef>) or an empty string.

=back

=head3 OUTPUT

The function will return an array reference containing hash references
with precursor information. Each precursor hash contain the following
key/value pairs:

    {
      precursor => "tca-...",       # precursor name
      len       => 121,             # length in basepairs
      seq       => "acug...",       # sequence of the precursor
      species   => "TCA",           # species tag
      matures   => [ ]              # matures belonging to the precursor
    }

The matures are represented as array containing hash references with
the following information:

    [
        {
            name  => "tca-...",      # mature name (non unique)
            start => 30,             # start coordinate in precursor sequence
            stop  => 50,             # end coordinate in precursor sequence
            seq   => "acug..."       # sequence
        },
        ...
    ]

=cut

sub parse_mirbase_dat
{
    my ($infile, $species) = @_;

    my @input = ();

    open(FH, "<", $infile) || die "Unable to open file '$infile' for reading: $!\n";
    while (<FH>)
    {
	# find a block starting with ^ID and ending with a line containing only //
	next unless (/^ID/);
	my $lines = $_;
	while (<FH>)
	{
	    $lines .= $_;
	    last if (/^\/\/$/);
	}

	my $block_content = _parse_mirbase_dat_block(\$lines, $species);
	push(@input, $block_content) if ($block_content);
    }
    close(FH) || die "Unable to close file '$infile' after reading: $!\n";

    # generate the mature sequences
    foreach my $entry (@input)
    {
	foreach my $mature (@{$entry->{matures}})
	{
	    $mature->{seq} = substr($entry->{seq}, $mature->{start}-1, $mature->{stop}-$mature->{start}+1);
	}
    }

    return \@input;
}

=pod

=head2 SUBROUTINE C<_parse_mirbase_dat_block()>

Internal function that will parse a single *.dat block.

=cut

sub _parse_mirbase_dat_block
{
    my ($ref_lines, $req_species) = @_;

    # first extract the ID line
    return unless ($$ref_lines =~ /^ID\s+(\S+)\s+(\S+);\s+(\S+);\s+(\S+);\s+(\d+)\s+BP.\s*$/m);
    my ($precursor, undef, undef, $species, $precursor_len) = ($1, $2, $3, $4, $5);

    return if ($req_species && uc($species) ne uc($req_species));

    # second we want to obtain the accession
    return unless ($$ref_lines =~ /^AC\s+(\S+);/m);
    my $accession = $1;

    # get the FT blocks for a mircoRNA
    my @matures = ();
    while ($$ref_lines =~ /^FT\s+miRNA\s+(\d+)\.\.(\d+).+?^FT\s+\/product="([^"]+)"/msg)
    {
	push(@matures, { start => $1, stop => $2, name => $3 });
    }

    # get the sequence
    return unless ($$ref_lines =~ /^SQ.+?$(.+)^\S/ms);
    my $seq = $1;
    $seq =~ s/\d+\s*$//mg;
    $seq =~ s/\s//g;

    return ( { precursor => $precursor, species => $species, len => $precursor_len, matures => \@matures, seq => $seq } );
}

=pod

=head2 SUBROUTINE C<parse_fasta()>

The subroutine will parse a fasta file and return a hash reference.

=head3 PARAMETER

=over 4

=item C<infile>

A filename of the input file

=back

=head3 OUTPUT

The function will return a hash reference containing the following
key/value pairs:

    {
      sequenceID => "actg...",
                                # sequenceID are the non-white characters behind the leading ">" in
                                # the fasta header. The corresponding value is the sequence with all
                                # whitespaces stripped out
    }

=cut

sub parse_fasta
{
    my ($infile) = @_;

    my %output;
    my $current_header;
    open(FH,"<",$infile) || die "Unable to open file '$infile' for reading: $!\n";
    while(<FH>)
    {
	chomp;
	if(/^>/)
	{
	    if ($_ =~ /^>(\S+)/)
	    {
		$current_header = $1;
	    } else {
		die "Something went wrong while parsing line '$_'\n";
	    }
	}
	else
	{
	    $output{$current_header} .= $_
	}
    }
    close(FH)|| die "Unable to close file '$infile' after reading: $!\n";
    return(\%output);
}

sub parse_mirdeep_novels
{
    my ($file, $cutoff)     = @_;
    my @result = @{_parse_mirdeep($file, $cutoff)->{novel}};

    @result = sort { $a->{digest} cmp $b->{digest} } (@result);

    return(\@result);
}

sub parse_mirdeep_known
{
    my ($file)     = @_;

    my @result = @{_parse_mirdeep($file, undef)->{known}};

    return(\@result);
}

sub _parse_mirdeep{
    my ($file, $cutoff)     = @_;
    my %result = (
	novel => [],
	known => [],
	);

    open(PM, "<", $file) || die "Unable to open file '$file': $!\n";
    while(<PM>){
	next unless (/^provisional/);

	my %seen = ();

	while(<PM>)
	{
	    chomp;

	    # leave the inner loop if the line is empty
	    last if (/^\s*$/);

	    my %dataset = ();

	    my @fieldnames = (
		"provisional_id",                         # provisional id
		"score",                                  # miRDeep2 score
		"probability",                            # estimated probability that the miRNA candidate is a true positive
		"rfam_alert",                             # rfam alert
		"total_read_count",                       # total read count
		"mature_count",                           # mature read count
		"loop_count",                             # loop read count
		"star_count",                             # star read count
		"randfold_significant",                   # significant randfold p-value
		"mirbase_mirna",                          # miRBase miRNA
		"mirbase_example",                        # example miRBase miRNA with the same seed
		"UCSC_browser",                           # UCSC browser
		"NCBI_blastn",                            # NCBI blastn
		"mature_seq",                             # consensus mature sequence
		"star_seq",                               # consensus star sequence
		"precursor_seq",                          # consensus precursor sequence
		"precursor_coordinates"                   # precursor coordinate
		);


	    @dataset{@fieldnames} = split("\t", $_);

	    next if ($cutoff && $dataset{score} < $cutoff);

	    # check if we need to switch mature and star
	    my $mature_idx	= index($dataset{precursor_seq},$dataset{mature_seq});
	    my $star_idx	= index($dataset{precursor_seq},$dataset{star_seq});

	    if ($star_idx == $mature_idx)
	    {
		die("mature and star sequence have the same position on hairpin\n"); # should never happen
	    }
	    # we need to switch if star is before mature
	    if($star_idx < $mature_idx){
		warn("Need to switch mature/star sequence for line '$_'\n");
		($dataset{mature_seq}, $dataset{star_seq}) = ($dataset{star_seq}, $dataset{mature_seq});
	    }

	    # calculate a checksum for "hairpin_seq|mature_seq|star_seq"
	    my $ctx = Digest::MD5->new;
	    my $current_seq = join("|", $dataset{precursor_seq}, $dataset{mature_seq}, $dataset{star_seq});
	    $ctx->add($current_seq);
	    my @checksum = unpack("S*", $ctx->digest); # digest is 16 Bytes... We are using the lowest 16 bit as unsigned number (ranging 0-65535) as digest
	    my $new_number = $checksum[-1];

	    $dataset{digest} = $new_number;
	    if (exists $seen{$new_number})
	    {
		# we found a collision
		# now check, if the sequences are identical
		my $former_hit = $seen{$new_number};
		my $former_seq = join("|", ($result{novel}[$former_hit]{precursor_seq}, $result{novel}[$former_hit]{mature_seq}, $result{novel}[$former_hit]{star_seq}));
		if ($former_seq eq $current_seq)
		{
		    # both sequences are identical, therefore we are assuming genomic copies
		    warn("Collision detected, but assuming to have found a genomic copy for line '$_'\n");
		} else {
		    die("Collision found, but sequences are different for line '$_'\n");
		}
	    } else {
		$seen{$new_number} = int(@{$result{novel}});
		push(@{$result{novel}}, \%dataset);
	    }
	}

	# identify the known block
	while(<PM>)
	{
	    last if (/^tag id/);
	}

	while(<PM>)
	{
	    chomp;

	    # leave the inner loop if the line is empty
	    last if (/^\s*$/);

	    my %dataset = ();

	    my @fieldnames = (
		"provisional_id",                         # provisional id
		"score",                                  # miRDeep2 score
		"probability",                            # estimated probability that the miRNA candidate is a true positive
		"rfam_alert",                             # rfam alert
		"total_read_count",                       # total read count
		"mature_count",                           # mature read count
		"loop_count",                             # loop read count
		"star_count",                             # star read count
		"randfold_significant",                   # significant randfold p-value
		"mirbase_mirna",                          # miRBase miRNA
		"mirbase_example",                        # example miRBase miRNA with the same seed
		"UCSC_browser",                           # UCSC browser
		"NCBI_blastn",                            # NCBI blastn
		"mature_seq",                             # consensus mature sequence
		"star_seq",                               # consensus star sequence
		"precursor_seq",                          # consensus precursor sequence
		"precursor_coordinates"                   # precursor coordinate
		);


	    @dataset{@fieldnames} = split("\t", $_);

	    # check if we need to switch mature and star
	    my $mature_idx	= index($dataset{precursor_seq},$dataset{mature_seq});
	    my $star_idx	= index($dataset{precursor_seq},$dataset{star_seq});

	    if ($star_idx == $mature_idx)
	    {
		die("mature and star sequence have the same position on hairpin\n"); # should never happen
	    }
	    # we need to switch if star is before mature
	    if($star_idx < $mature_idx){
		warn("Need to switch mature/star sequence for line '$_'\n");
		($dataset{mature_seq}, $dataset{star_seq}) = ($dataset{star_seq}, $dataset{mature_seq});
	    }

	    push(@{$result{known}}, \%dataset);
	}

	# we will arrive here, if we already parsed the novel and
	# known block and can therefore leave the outer loop as
	# well
	last;
    }

    close(PM) || die "Unable to close file '$file': $!\n";

    return(\%result);
}


1;
