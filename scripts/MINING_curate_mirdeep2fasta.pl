#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

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
die "Need to specify --cutoff cutoffscore\n" unless (defined $cutoff);
die "Need to specify --matureout file\n" unless (defined $mature_file);
die "Need to specify --hairpinout file\n" unless (defined $hairpin_file);
die "Need to specify --species three-letter-species-code\n" unless (defined $species);

# get novel miRNAs above threshold
my %novel_hash	= %{&parse_mirdeep($csv_file,$cutoff)};

open(MATURE,">",$mature_file)   || die "Unable to open file '$mature_file': $!\n";
open(HAIRPIN,">",$hairpin_file) || die "Unable to open file '$hairpin_file': $!\n";

foreach(keys %novel_hash){
	my $novel_count	= $_;
	my @novel_split	= @{$novel_hash{$novel_count}};
	my $mature	= uc($novel_split[7]);
	my $star	= uc($novel_split[8]);
	my $hairpin	= uc($novel_split[9]);
	
	my $mature5p;
	my $mature3p;
	
	my $mature_idx	= index($hairpin,$mature);
	my $star_idx	= index($hairpin,$star);
	
	if($mature_idx < $star_idx){
		$mature5p=$mature;
		$mature3p=$star;
	}
	elsif($star_idx < $mature_idx){
		$mature5p=$star;
		$mature3p=$mature;
	}
	else{
		die("mature and star sequence have the same position on hairpin\n"); # should never happen
	}	

	$mature5p	=~ s/U/T/g;
	$mature3p	=~ s/U/T/g;
	$hairpin	=~ s/U/T/g;

	my $header	= sprintf(">%s-new-%d", $species, $novel_count);

	print HAIRPIN "$header\n$hairpin\n";
	print MATURE "$header-5p\n$mature5p\n$header-3p\n$mature3p\n";
}
#close(CSV) || die;

close(MATURE) || die "Unable to close file '$mature_file': $!\n";
close(HAIRPIN)|| die "Unable to close file '$hairpin_file': $!\n";

###########################################################################

# skip all parts but novel section
sub parse_mirdeep{
        my $pm_file     = $_[0];
        my $pm_cutoff   = $_[1];
        my %pm_hash;
        my $pm_bool     = 0;
	my $pm_count	= 1;
        open(PM,"<",$pm_file) || die;
        while(<PM>){
                chomp;
                my $pm_line    = $_;
                if (/^$/){
                        $pm_bool = 0;
                        next;
                }
                if(/^provisional/){
                        $pm_bool = 1;
                        next;
                }
                if($pm_bool == 1){
                        my @pm_split            = split("\t",$pm_line);
			my ($pm_prov_id,$pm_score,undef,undef,undef,$pm_mature_count,$pm_loop_count,$pm_star_count,undef,undef,$pm_mirbase_ref,undef,undef,$pm_mature_seq,$pm_star_seq,$pm_precursor_seq,$pm_precursor_coord) = split("\t",$pm_line);
                        next if ($pm_score < $pm_cutoff);
                        my @pm_list     = ($pm_file,$pm_prov_id,$pm_score,$pm_mature_count,$pm_loop_count,$pm_star_count,$pm_mirbase_ref,$pm_mature_seq,$pm_star_seq,$pm_precursor_seq,$pm_precursor_coord);
                        $pm_hash{$pm_count}=\@pm_list;
			$pm_count++;
                }
        }
        close(PM) || die;
        return(\%pm_hash);
}

