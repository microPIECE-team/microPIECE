#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Digest::MD5;

my $csv_file;
my $cutoff;
my $mature_file;
my $hairpin_file;
my $species;

GetOptions(
    "csv=s"        => \$csv_file,
    "cutoff=i"     => \$cutoff,
    "matureout=s"  => \$mature_file,
    "hairpinout=s" => \$hairpin_file,
    "species=s"    => \$species,
) || die;

die "Need to specify --csv file\n" unless (defined $csv_file);
die "Input file does not exist\n" unless (-e $csv_file);
die "Need to specify --cutoff cutoffscore\n" unless (defined $cutoff);
die "Need to specify --matureout file\n" unless (defined $mature_file);
die "Mature file exists and will not be overwritten!\n" if (-e $mature_file);
die "Need to specify --hairpinout file\n" unless (defined $hairpin_file);
die "Hairpin file exists and will not be overwritten!\n" if (-e $hairpin_file);
die "Need to specify --species three-letter-species-code\n" unless (defined $species);

# get novel miRNAs above threshold
my $novels = parse_mirdeep($csv_file, $cutoff);

open(MATURE,">",$mature_file)   || die "Unable to open file '$mature_file': $!\n";
open(HAIRPIN,">",$hairpin_file) || die "Unable to open file '$hairpin_file': $!\n";

for(my $novel_count = 1; $novel_count <= @{$novels}; $novel_count++)
{
    my $novel = $novels->[$novel_count-1];
    
    my $mature5p = uc($novel->{mature_seq});
    my $mature3p = uc($novel->{star_seq});
    my $hairpin	 = uc($novel->{precursor_seq});

    $mature5p	=~ s/U/T/g;
    $mature3p	=~ s/U/T/g;
    $hairpin	=~ s/U/T/g;

    my $header	= sprintf(">%s-new-%d", $species, $novel->{digest});

    print HAIRPIN $header, "\n", $hairpin, "\n";
    print MATURE  $header, "-5p\n", $mature5p, "\n", $header, "-3p\n", $mature3p, "\n";
}
close(MATURE) || die "Unable to close file '$mature_file': $!\n";
close(HAIRPIN)|| die "Unable to close file '$hairpin_file': $!\n";

###########################################################################

# skip all parts but novel section
sub parse_mirdeep{
        my ($file, $cutoff)     = @_;
        my @result = ();
	my %seen   = ();

        open(PM, "<", $file) || die "Unable to open file '$file': $!\n";
        while(<PM>){
	    next unless (/^provisional/);

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

		next if ($dataset{score} < $cutoff);

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
		    my $former_seq = join("|", ($result[$former_hit]{precursor_seq}, $result[$former_hit]{mature_seq}, $result[$former_hit]{star_seq}));
		    if ($former_seq eq $current_seq)
		    {
			# both sequences are identical, therefore we are assuming genomic copies
			warn("Collision detected, but assuming to have found a genomic copy for line '$_'\n");
		    } else {
			die("Collision found, but sequences are different for line '$_'\n");
		    }
		} else {
		    $seen{$new_number} = int(@result);
		    push(@result, \%dataset);
		}
	    }

	    # we will arrive here, if we already parsed the novel
	    # block and can therefore leave the outer loop as well
	    last;
	}

	close(PM) || die "Unable to close file '$file': $!\n";

	@result = sort { $a->{digest} cmp $b->{digest} } (@result);

        return(\@result);
}
