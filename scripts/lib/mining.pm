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

sub export_fasta
{
    my ($file_mature, $file_hairpin, $data) = @_;

    # use a dummy, if a file is not specified
    unless (defined $file_mature)
    {
	$file_mature = \do{my $mature = ""};
    }
    unless (defined $file_hairpin)
    {
	$file_hairpin = \do{my $hairpin = ""};
    }

    open(MATURE, ">", $file_mature) || die "Unable to open file '$file_mature': $!\n";
    open(HAIRPIN, ">", $file_hairpin) || die "Unable to open file '$file_hairpin': $!\n";

    my %seen = ( matures => {}, precursor => {} );

    foreach my $entry (@{$data})
    {
	# print the precursor
	my $seq = $entry->{seq};
	$seq = uc($seq);    # upper case
	$seq =~ tr/U/T/;    # replace Us by Ts
	my $fasta_block = ">".$entry->{precursor}."\n".$seq."\n";
	unless (exists $seen{precursor}{$fasta_block})
	{
	    print HAIRPIN $fasta_block;
	}
	$seen{precursor}{$fasta_block}++;

	# print the matures
	foreach my $mature (@{$entry->{matures}})
	{
	    # start and stop have to exist, seq not, therefore, we
	    # will generate the sequence on the fly
	    my $seq = substr($entry->{seq}, $mature->{start}-1, $mature->{stop}-$mature->{start}+1);
	    $seq = uc($seq);    # upper case
	    $seq =~ tr/U/T/;    # replace Us by Ts
	    my $fasta_block = ">".$mature->{name}."\n".$seq."\n";
	    unless (exists $seen{mature}{$fasta_block})
	    {
		print MATURE $fasta_block;
	    }
	    $seen{mature}{$fasta_block}++;
	}
    }

    close(HAIRPIN) || die "Unable to close file '$file_hairpin': $!\n";
    close(MATURE) || die "Unable to close file '$file_mature': $!\n";
}

sub export_mirbase_data
{
    my ($file, $data) = @_;

    open(FH, ">>", $file) || die "Unable to open file '$file': $!\n";
    foreach my $entry (@{$data})
    {
	my $buffer = "";

	$buffer .= sprintf("ID   %s    standard; RNA; %s; %d BP.\nXX\n", $entry->{precursor}, $entry->{species}, length($entry->{seq}));
	if (exists $entry->{description})
	{
	    $buffer .= sprintf("DE   %s\nXX\n", $entry->{description});
	}
	$buffer .= "FH   Key             Location/Qualifiers\nFH\n";
	foreach my $mature (@{$entry->{matures}})
	{
	    $buffer .= sprintf("FT   miRNA           %d..%d\n", $mature->{start}, $mature->{stop});
	    $buffer .= sprintf('FT                   /product="%s"'."\n", $mature->{name});
	}
	$buffer.= "XX\n";

	# add sequence information
	my $counterA = $entry->{seq} =~ tr/Aa/Aa/;
	my $counterC = $entry->{seq} =~ tr/Cc/Cc/;
	my $counterG = $entry->{seq} =~ tr/Gg/Gg/;
	my $counterT = $entry->{seq} =~ tr/Tt/Tt/;
	my $counterOthers = length($entry->{seq})-$counterA-$counterC-$counterG-$counterT;

	$buffer .= sprintf("SQ   Sequence %d BP; %d A; %d C; %d G; %d T; %d other;\n", length($entry->{seq}), $counterA, $counterC, $counterG, $counterT, $counterOthers);
	for(my $i=0; $i<length($entry->{seq}); $i+=60)
	{
	    my $subseq = substr($entry->{seq}, $i, 60);
	    my $end    = $i+length($subseq);
	    # group of 10 characters
	    $subseq =~s/(.{10})/$1 /g;
	    $subseq = "     ".$subseq;
	    my $numberlength=80-length($subseq);
	    $buffer .= sprintf("%s% ".$numberlength."d\n", $subseq, $end);
	}
	$buffer .= "//\n";

	print FH $buffer;
    }
    close(FH) || die "Unable to close file '$file': $!\n";
}

sub fix_hairpin
{
    my ($structure, $sequence) = @_;

    # split the structure into its lines
    my @lines = split(/\n/, $structure);

    # and each line into its characters
    @lines = map { [split("", $_)] } (@lines);

    # split the input sequence into its characters
    my @seq = split("", $sequence);

    my $width = int(@{$lines[0]});

    # now go through the lines 1 and 2 and substitute each character
    # found, with the next character from our input sequence and
    # second test lines 4 and 5 for the same condition and use the
    # last character from the sequence as substitution character
    for(my $x=0; $x<$width; $x++)
    {
	# test if y=0 or y=1 contains the character
	my $y;
	if ($lines[0][$x] !~ /[-\s]/)
	{
	    $y = 0;
	} elsif ($lines[1][$x] !~ /[-\s]/)
	{
	    $y = 1;
	}

	# set the next character from our input sequence
	$lines[$y][$x] = shift @seq if (defined $y);

	# test if y=3 or y=4 contains the character
	$y=undef;
	if ($lines[3][$x] !~ /[-\s]/)
	{
	    $y = 3;
	} elsif ($lines[4][$x] !~ /[-\s]/)
	{
	    $y = 4;
	}

	# set the next character from our input sequence
	$lines[$y][$x] = pop @seq if (defined $y)
    }

    # if sequence still contain a letter, it should be on line 3 at
    # the very last position
    if (@seq)
    {
	if ($lines[2][$width-1] !~ /[-\s]/)
	{
	    $lines[2][$width-1] = shift @seq;
	}
    }

    # last join each line and all lines together
    return join("\n", map { join("", @{$_}) } (@lines));
}

1;
