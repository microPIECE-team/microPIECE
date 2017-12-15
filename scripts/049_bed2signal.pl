#!/usr/bin/perl
use strict;
use warnings;
# takes ff bedfile and extracts signal positions
# signal strength defined by user


my $bed_file		= $ARGV[0];
my $signal_strength	= $ARGV[1];

open(BED,"<",$bed_file) || die;
while(<BED>){
	chomp;

	# ignore comment lines and pipe them through
	if($_ =~ /^#/)
	{
	    print "$_\n";
	    next;
	}

	# parse the bed line
	my $bed_line	= $_;	# chr	start	stop	info	.	strand
	my ($bed_chr, $bed_start, $bed_stop, $bed_info, undef, $bed_strand)	= split("\t",$bed_line);
	my %data_from_bed_info = ();
	foreach my $key_value (split(";",$bed_info))           # length=XY;counts=pos_1sum/pos1_rep1/pos2_rep2/.../pos1_repN/,pos2_sum/pos2_rep1/...
	{
	    my ($key, $value) = split("=", $key_value,2);
	    $data_from_bed_info{$key} = $value;
	}

	# split the count data
	$data_from_bed_info{counts} = [ map { [ split('/', $_) ] } (split(",", $data_from_bed_info{counts})) ];

	# go through the bed counts and keep stretches which are
	# supported by $signal_strength or more peaks and store the
	# coordinates in new array @subregions
	my @subregions = ();

	for(my $i=0; $i<@{$data_from_bed_info{counts}}; $i++)
	{
	    if ($data_from_bed_info{counts}[$i][0] >= $signal_strength)
	    {
		my $start = $i;
		my $stop  = $i;
		for(my $j=$i+1; $j<@{$data_from_bed_info{counts}}; $j++)
		{
		    if ($data_from_bed_info{counts}[$j][0] < $signal_strength)
		    {
			last;
		    } else {
			$stop = $j;
		    }
		}

		$i = $stop; # ensure we are starting after the region

		push(@subregions, { start => $start, stop => $stop } );
	    }
	}

	# print subregions as new bed lines
	foreach my $subregion (@subregions)
	{
	    my $new_bed_chr    = $bed_chr;                       # should be the same
	    my $new_bed_start  = $bed_start+$subregion->{start}; # shift the new start
	    my $new_bed_stop   = $bed_start+$subregion->{stop};  # shift the new stop
	    my $new_bed_info   = get_new_info_string(\%data_from_bed_info, $subregion);
	    my $new_bed_strand = $bed_strand;                    # should be the same

	    print join("\t",
		       $new_bed_chr,
		       $new_bed_start,
		       $new_bed_stop,
		       $new_bed_info,
		       ".",
		       $new_bed_strand
		), "\n";

	}
}
close(BED) || die;


sub get_new_info_string
{
    my ($ref_data, $ref_subregion) = @_;

    my @output = ();

    # create a string containing all information, but counts and length
    foreach my $key (keys %{$ref_data})
    {
	next if ($key eq "counts" || $key eq "length");

	push(@output, $key."=".$ref_data->{$key});
    }

    # get the counts for the subregion
    my @new_counts = map { $ref_data->{counts}[$_] } ($ref_subregion->{start}..$ref_subregion->{stop});
    push(@output, "counts=".join(",", map { join("/", @{$_}) } (@new_counts)));

    # get the new length
    push(@output, "length=".int(@new_counts));

    # add a new key-value-pair to indicate that string as substring
    push(@output, "subregion=".join(",", ($ref_subregion->{start},$ref_subregion->{stop})));

    # add a new key-value-pair to indicate the original length
    push(@output, "originallength=".$ref_data->{length});

    return join(";", @output);
}
