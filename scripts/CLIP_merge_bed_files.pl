#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;

=pod

=head1 NAME

046_merge_bed_files.pl - merging bed files, while maintaining count information

=head1 DESCRIPTION

This script is used to merge bed files into a single bed file, while
keeping the information about the input counts in column #4 of the
resulting bed file.

=head1 SYNOPSIS

   046_merge_bed_files.pl [options] --input inputclass=filename1,filename2

   Options:
     --help
     --output
     --overwrite
     --log

   # Example defining to classes with two files each and output to merged.bed
   046_merge_bed_files.pl
      --input 24h=file1.bed,file2.bed \
      --input 72h=file3.bed,file4.bed \
      --output merged.bed

=head1 OPTIONS

=over 8

=item C<--input>

Specifies the input classes. For a single class, at least one file
need to be specified. Multiple files are supported. In that case use a
comma as delimiter between file names.

Counts for each input class are maintained.

Specifying at least one input class is mandatory.

=item C<--output>

Specifies the output file. Will be a bed file, with columns 1-3
default bed format, followed by a count/length field, an undefiled
column (.), and the strand field.

=item C<--overwrite>

By default, an existing output file will not overwritten, but will
kill the script. Using --overwrite will force to overwrite existing
output files.

=item C<--log>

Specifies the location of the log file. If no log file is set, logging
information are written to STDERR, if a file is set, logging to
STDERR/that file will performed.

=item C<--help>

Print that help message

=back

=cut

my %input     = ();
my $output    = "-"; # default is standard out
my $log       = "STDERR"; # default is standard error, if a value is specified, the messages will be printed to STDERR and that file
my $overwrite = 0;
my $help      = 0;

use Getopt::Long;
use Log::Log4perl qw(:easy);

GetOptions(
    "input=s%"  => \%input,
    "output=s"  => \$output,
    "log=s"     => \$log,
    "overwrite" => \$overwrite,
    "help"       => \$help
    ) || pod2usage(2);

# check that at least one class was specified via input parameter
unless (keys %input >= 1)
{
    pod2usage(2);
}

pod2usage(1) if $help;

# split input keys into seperate files
foreach my $key (keys %input)
{
    $input{$key} = [ split(/,/, $input{$key}) ];
}

my @log4perl_init = (
    {
	level  => $DEBUG,
	file   => "STDERR",
	layout => '[%d] (%p) %m%n'
    }
    );

if ($log ne "STDERR")
{
    push(@log4perl_init,
	 {
	     level  => $DEBUG,
	     file   => ">>$log",
	     layout => '[%d] (%p) %m%n'
	 }
	);
}

Log::Log4perl->easy_init(@log4perl_init);

# check if output is redirected to a file
my $fh = *STDOUT;  # default output should go to STDOUT
if ($output ne "-")
{
    # test if output file exists
    if (-e $output && ! $overwrite)
    {
	LOGDIE("Output file '$output' exists!");
    } else {
	open($fh, ">", $output) || LOGDIE("Unable to open output file '$output': $!");
    }
}

my @genome = ();
my %chromosomes = ();
my %strands    = ();
my %conditions = ();


foreach my $key (keys %input)
{
    foreach my $file (@{$input{$key}})
    {
	WARN("Working on file '$file'");
	open(FH, "<", $file) || LOGDIE("Unable to open file '$file': $!");
	while(<FH>)
	{
	    # go through the bed files
	    chomp;

	    # NW_001809801.1	79510	79534	XM_001647792.1	1	+
	    # NW_001809801.1	79527	79534	XM_001647792.1	1	+
	    # NW_001809801.1	248783	248807	XM_001647796.1	1	-
            # NW_001809801.1	533942	533950	XM_001647802.1	1	-
            # NW_001809801.1	2779750	2779790	XM_001647839.1	1	+

	    my ($chromosome, $start, $stop, $mrna_ids, undef, $strand) = split(/\t/, $_);

	    # check if chromosome is already known
	    unless (exists $chromosomes{$chromosome})
	    {
		$chromosomes{$chromosome} = int(keys %chromosomes);
	    }
	    $chromosome = $chromosomes{$chromosome};

	    # check if strand is already known
	    unless (exists $strands{$strand})
	    {
		$strands{$strand} = int(keys %strands);
	    }
	    $strand = $strands{$strand};

	    # check if condition/key is already known
	    unless (exists $conditions{$key})
	    {
		$conditions{$key} = int(keys %conditions);
	    }
	    my $condition = $conditions{$key};

	    # for each chromosomal position we need to consider:
	    # Counter + strand
	    # Counter - strand

	    # due to bed stop field is 0-based, but exclusive, we
	    # need to continue as long as $i is less then stop-field
	    for (my $i=$start; $i<$stop; $i++)
	    {
		$genome[$chromosome][$strand][$i][$condition]++;
	    }
	}
	close(FH) || LOGDIE("Unable to close file '$file': $!");
    }
}

# print the output
my @conditions_ordered = sort (keys %conditions);
print $fh "# Conditional counts are printed in the following order: ", join(", ", ("total", @conditions_ordered)), "\n";

foreach my $chromosome_key (sort keys %chromosomes)
{
    my $chromosome = $chromosomes{$chromosome_key};
    foreach my $strand_key (sort keys %strands)
    {
	my $strand = $strands{$strand_key};

	# check if the strand and the chromosome are existing
	unless (defined $genome[$chromosome] && ref($genome[$chromosome]) eq "ARRAY")
	{
	    # chromosomes should be always defined!
	    LOGDIE("Missing entry for chromosome '$chromosome_key'");
	}
	unless (defined $genome[$chromosome][$strand] && ref($genome[$chromosome][$strand]) eq "ARRAY")
	{
	    # if no feature is annotated on the strand if could be missing
	    WARN("No feature on chromosome '$chromosome_key' for ($strand_key)-strand.");
	    next;
	}

	my $start = -1;
	for (my $i=0; $i<@{$genome[$chromosome][$strand]}; $i++)
	{
	    # check if we can skip this position
	    unless (defined $genome[$chromosome][$strand][$i] && ref($genome[$chromosome][$strand][$i]) eq "ARRAY")
	    {
		next;
	    }

	    my $counts_on_position = get_counts_for_position(\@genome, $chromosome, $strand, $i, \%conditions);

	    if ($counts_on_position->{total} > 0)
	    {
		# we found a new block
		$start = $i;
		my $stop = -1;

		my @counts = ($counts_on_position);

		for (my $j=$i+1; $j<@{$genome[$chromosome][$strand]}; $j++)
		{
		    my $counts_on_inner_position = get_counts_for_position(\@genome, $chromosome, $strand, $j, \%conditions);

		    if ($counts_on_inner_position->{total} > 0)
		    {
			$stop = $j;
			push(@counts, $counts_on_inner_position);
		    } else {
			# due to bed stop field is 0-based, but exclusive, we
			# need to increase stop coordinate by 1
			$stop++;
			last;
		    }
		}
		if ($stop == -1)
		{
		    # due to bed stop field is 0-based, but exclusive, we
		    # need to use the array size
		    $stop = int(@{$genome[$chromosome][$strand]});
		}

		print $fh join("\t", ($chromosome_key, $start, $stop, sprintf("length=%d;counts=%s", @counts+0, generate_cigar_like_string(\@counts, ["total", @conditions_ordered])), ".", $strand_key)), "\n";
		$i = $stop;
		$start = -1; $stop = -1;
	    }
	}
	$genome[$chromosome][$strand] = undef; # reduce memory footprint
    }
    $genome[$chromosome] = undef; # reduce memory footprint
}

sub generate_cigar_like_string
{
    my ($ref_counts, $ref_ordered_conditions) = @_;

    my @output = ();

    foreach my $counts (@{$ref_counts})
    {
	push(@output, join("/", map {$counts->{$_}} (@{$ref_ordered_conditions})));
    }

    return join(",", @output);
}

sub get_counts_for_position
{
    my ($ref_genome, $chr, $str, $pos, $ref_conditions) = @_;

    my %counts = ( total => 0 );

    foreach my $condition (keys %{$ref_conditions})
    {
	my $val = $ref_genome->[$chr][$str][$pos][$ref_conditions->{$condition}];
	$counts{$condition} = (defined $val) ? $val : 0;
	$counts{total} += $counts{$condition};
    }

    return \%counts;
}
