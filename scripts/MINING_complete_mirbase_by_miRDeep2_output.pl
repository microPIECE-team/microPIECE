#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/lib";

use mining;

my $mirdeep_csv;
my $mirbase_dat;
my $species = "";

GetOptions(
	"mirdeep_out=s"		=>\$mirdeep_csv,
	"mirbase_dat=s"	        =>\$mirbase_dat,
        "species=s"             =>\$species) || die;

my $mirbase_dat_content = mining::parse_mirbase_dat($mirbase_dat, $species);

my %mature2precursorIndex = ();
for(my $i=0; $i<@{$mirbase_dat_content}; $i++)
{
    foreach my $mature (@{$mirbase_dat_content->[$i]{matures}})
    {
	push(@{$mature2precursorIndex{$mature->{name}}}, $i);
    }
}

my $mirdeep_content = mining::parse_mirdeep_known($mirdeep_csv);

foreach my $mirdeep_line (@{$mirdeep_content})
{
    unless (exists $mature2precursorIndex{$mirdeep_line->{mirbase_mirna}})
    {
	die "Unable to find an entry for $mirdeep_line->{mirbase_mirna}\n";
    }

    # get a list of potential precursor
    my @precursor = @{$mature2precursorIndex{$mirdeep_line->{mirbase_mirna}}};

    # filter for precursors without two mature sequence
    my @precursor_2_complete = grep { int(@{$mirbase_dat_content->[$_]{matures}}) < 2 } (@precursor);

    # for the precursors_2_complete, we want add another mature from mirdeep
    if (@precursor_2_complete)
    {

	# check for a minimum of 10 reads each for mature/star sequence or sum of both > 0 and both minimum 5
	if (
	    ($mirdeep_line->{mature_count} >= 10 && $mirdeep_line->{star_count} >= 10)
	    ||
	    ($mirdeep_line->{mature_count}+$mirdeep_line->{star_count}>=100 && $mirdeep_line->{mature_count}>=5 && $mirdeep_line->{star_count}>=5)
	    )
	{
	    warn(sprintf("Mature/star sequences for '%s' do have enough counts (%d/%d) to fulfil requirement for adding mature sequence\n", $mirdeep_line->{mirbase_mirna}, $mirdeep_line->{mature_count}, $mirdeep_line->{star_count}));
	} else {
	    warn(sprintf("Mature/star sequences for '%s' do not have enough counts (%d/%d) to fulfil requirement for adding mature sequence\n", $mirdeep_line->{mirbase_mirna}, $mirdeep_line->{mature_count}, $mirdeep_line->{star_count}));
	    next;
	}

	warn("Adding missing mature sequences for mirDeep2 mature ".$mirdeep_line->{mirbase_mirna}."\n");

	if (@precursor > 1 && @precursor_2_complete != @precursor)
	{
	    warn("Multiple potential precursors with some not having 2 matures for '".$mirdeep_line->{mirbase_mirna}."\n");
	}

	foreach my $precursor (@precursor_2_complete)
	{
	    my $entry = $mirbase_dat_content->[$precursor];

	    my $matures = $entry->{matures};

	    # add the missing mature
	    my $new_mature_entry = {
		name  => $entry->{precursor},
		start => index($entry->{seq}, $mirdeep_line->{star_seq}),
		stop  => index($entry->{seq}, $mirdeep_line->{star_seq})+length($mirdeep_line->{star_seq})-1,
		seq   => $mirdeep_line->{star_seq}
	    };

	    # check the positions:
	    if ($matures->[0]{name} =~ /-3p$/)
	    {
		if ($new_mature_entry->{start} < $matures->[0]{start})
		{
		    # as expected, the new gets -5p
		    $new_mature_entry->{name} .= "-5p";
		} else {
		    die "Unexpected new mature in line '".join(",", $mirdeep_line)."\n";
		}
	    } elsif ($matures->[0]{name} =~ /-5p$/) {
		if ($new_mature_entry->{start} > $matures->[0]{start})
		{
		    # as expected, the new gets -3p
		    $new_mature_entry->{name} .= "-3p";
		} else {
		    die "Unexpected new mature in line '".join(",", $mirdeep_line)."\n";
		}
	    } else {
		warn "Unable to determine if 3p/5p for name ".$matures->[0]{name}.", but will rename both matures correctly\n";
		if ($new_mature_entry->{start} > $matures->[0]{start})
		{
		    $new_mature_entry->{name} .= "-3p";
		    $matures->[0]{name}       .= "-5p";
		} else {
		    $new_mature_entry->{name} .= "-5p";
		    $matures->[0]{name}       .= "-3p";
		}
	    }

	    push(@{$matures}, $new_mature_entry);

	    # add to the description to the entry
	    $entry->{description} .= " Added mature sequence fulfilling the mirbase high confidence criteria";
	}
    }
}

my $output = "";
mining::export_mirbase_dat(\$output, $mirbase_dat_content);

print $output;
