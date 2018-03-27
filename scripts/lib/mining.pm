package mining;

use strict;
use warnings;

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

	my $block_content = parse_mirbase_dat_block(\$lines, $species);
	push(@input, $block_content) if ($block_content);
    }
    close(FH) || die "Unable to close file '$infile' after reading: $!\n";

    return \@input;
}

sub parse_mirbase_dat_block
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
