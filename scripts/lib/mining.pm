package mining;

use strict;
use warnings;

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


1;
