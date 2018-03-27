#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Digest::MD5;

use FindBin;
use lib "$FindBin::Bin/lib";

use mining;

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
my $novels = mining::parse_mirdeep_novels($csv_file, $cutoff);

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
