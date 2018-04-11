#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;

=pod

=head1 NAME

CLIP_merge_bed_files.pl - merging bed files, while maintaining count information

=head1 DESCRIPTION

This script is used to merge bed files into a single bed file, while
keeping the information about the input counts in column #4 of the
resulting bed file.

=head1 SYNOPSIS

   CLIP_merge_bed_files.pl [options] --input inputclass=filename1,filename2

   Options:
     --help
     --output
     --overwrite
     --log

   # Example defining to classes with two files each and output to merged.bed
   CLIP_merge_bed_files.pl
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

my %chromosomes = ();
my %strands    = ();

my %input_sorted = ();
my $max_length = 0;

foreach my $key (keys %input)
{
    foreach my $file (@{$input{$key}})
    {
	WARN("Working on file '$file'");
	WARN("Sorting file '$file'");
	open(FH, "<", $file) || LOGDIE("Unable to open file '$file': $!");
	my @temp = <FH>;
	@temp = map { [ split("\t", $_) ] } (@temp);

	# sort by chromosome, followed by strand, followed by start coordinate, followed by stop coordinate
	@temp = sort { $a->[0] cmp $b->[0] || $a->[5] cmp $b->[5] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } (@temp);

	# estimate chromosomes and their maximum coordinate
	foreach my $set (@temp)
	{
	    $chromosomes{$set->[0]} = $set->[2] unless (exists $chromosomes{$set->[0]} && $chromosomes{$set->[0]} >= $set->[2]);
	    $max_length = $set->[2] if ($max_length < $set->[2]);

	    # store strands
	    $strands{$set->[5]}++;
	}

	push(@{$input_sorted{$key}{dat}} , \@temp);
	push(@{$input_sorted{$key}{pos}} , 0);
    }
}

my @chromosome_order = sort {$a cmp $b} (keys %chromosomes);
WARN(sprintf("Found a total of %d chromosomes with a maximum length (covered) of %d bp. Expected chromosome order: %s", int(@chromosome_order), $max_length, join(", ", map {"'$_'"} (@chromosome_order))));

# map condition to number
my @conditions_ordered = sort (keys %input_sorted);
my %condition2number;
@condition2number{@conditions_ordered} = (0..(int(@conditions_ordered)-1));

# map strand to number
my @strands_ordered = sort (keys %strands);
my %strand2number;
@strand2number{@strands_ordered} = (0..(int(@strands_ordered)-1));

print $fh "# Conditional counts are printed in the following order: ", join(", ", ("total", @conditions_ordered)), "\n";

foreach my $chromosome (@chromosome_order)
{

    my @genome = ();

    WARN(sprintf("Working on chromosome %s", $chromosome));

    foreach my $condition (keys %input_sorted)
    {
	for(my $entry=0; $entry <@{$input_sorted{$condition}{dat}}; $entry++)
	{
	    next if ($input_sorted{$condition}{pos}[$entry] >= @{$input_sorted{$condition}{dat}[$entry]});

	    for (my $i=$input_sorted{$condition}{pos}[$entry]; $i<@{$input_sorted{$condition}{dat}[$entry]}; $i++)
	    {
		# check if the chromosome is still correct
		my ($chr, $start, $stop, $mrna_ids, $score, $strand) = @{$input_sorted{$condition}{dat}[$entry][$i]};  # a single line from a input bed
		# NW_001809801.1	79510	79534	XM_001647792.1	1	+
		# NW_001809801.1	79527	79534	XM_001647792.1	1	+
		# NW_001809801.1	248783	248807	XM_001647796.1	1	-
		# NW_001809801.1	533942	533950	XM_001647802.1	1	-
		# NW_001809801.1	2779750	2779790	XM_001647839.1	1	+

		# if the chromosome is not expected, we will store the position and go on
		if ($chr ne $chromosome)
		{
		    $input_sorted{$condition}{pos}[$entry] = $i;
		    last;
		}

		# for each chromosomal position we need to consider:
		# Counter + strand
		# Counter - strand

		# due to bed stop field is 0-based, but exclusive, we
		# need to continue as long as $i is less then stop-field
		for (my $j=$start; $j<$stop; $j++)
		{
		    $genome[$strand2number{$strand}][$j][$condition2number{$condition}]++;
		}
	    }
	}
    }

    foreach my $strand (@strands_ordered)
    {
	unless (defined $genome[$strand2number{$strand}] && ref($genome[$strand2number{$strand}]) eq "ARRAY")
	{
	    # if no feature is annotated on the strand if could be missing
	    WARN("No feature on chromosome '$chromosome' for ($strand)-strand.");
	    next;
	}

	my $start = -1;
	for (my $i=0; $i<@{$genome[$strand2number{$strand}]}; $i++)
	{
	    # check if we can skip this position
	    unless (defined $genome[$strand2number{$strand}][$i] && ref($genome[$strand2number{$strand}][$i]) eq "ARRAY")
	    {
		next;
	    }

	    my $counts_on_position = get_counts_for_position($genome[$strand2number{$strand}][$i], \@conditions_ordered);

	    if ($counts_on_position->{total} > 0)
	    {
		# we found a new block
		$start = $i;
		my $stop = $i;

		my @counts = ($counts_on_position);

		for (my $j=$i+1; $j<@{$genome[$strand2number{$strand}]}; $j++)
		{
		    my $counts_on_inner_position = get_counts_for_position($genome[$strand2number{$strand}][$j], \@conditions_ordered);

		    if ($counts_on_inner_position->{total} > 0)
		    {
			$stop = $j;
			push(@counts, $counts_on_inner_position);
		    } else {
			last;
		    }
		}

		# due to bed stop field is 0-based, but exclusive, we
                # need to increase stop coordinate by 1
		print $fh join("\t", (
				   $chromosome,
				   $start,
				   $stop+1,
				   sprintf("length=%d;counts=%s",
					   @counts+0,
					   generate_cigar_like_string(\@counts, ["total", @conditions_ordered]
					   )
				   ),
				   ".",
				   $strand
			       )
		    ), "\n";
		$i = $stop; ## will be increased by the for loop, therefore next iteration starts at $stop+1
	    }
	}

    }
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
    my ($ref_pos, $ref_conditions_ordered) = @_;

    my %counts = ( total => 0 );

    for(my $i=0; $i<@{$ref_conditions_ordered}; $i++)
    {
	my $val       = $ref_pos->[$i];
	my $condition = $ref_conditions_ordered->[$i];

	$counts{$condition}  = (defined $val) ? $val : 0;
	$counts{total}      += $counts{$condition};
    }

    return \%counts;
}
