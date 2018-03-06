#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $precursor_file;
my $mature_file;
my $species;		# tca
my $organism_file;	# organisms.txt from mirbase
my $out_species_mature;
my $out_species_precursor;
my $out_nonspecies_mature;

GetOptions(
	"organisms=s"		=> \$organism_file,
	"species=s"		=> \$species,
        "precursor_file=s"      => \$precursor_file,
        "mature_file=s" 	=> \$mature_file,
	"outmature=s"		=> \$out_species_mature,
	"outprecursor=s"	=> \$out_species_precursor,
	"outnonspeciesmature=s"	=> \$out_nonspecies_mature,
) || die;

die "Need to specify precursor_file via --precursor_file *filename* option\n" unless (defined $precursor_file && -e $precursor_file);
die "Need to specify mature_file via --mature_file *filename* option\n" unless (defined $mature_file && -e $mature_file);
die "Need to specify organism_file via --organism *filename* option\n" unless (defined $organism_file && -e $organism_file);
die "Need to specify the species via --species *speciestag* option\n" unless (defined $species);
die "Need to specify a non existing mature output file via --outmature *filename* option\n" if ((! defined $out_species_mature) || (-e $out_species_mature));
die "Need to specify a non existing precursor output file via --outprecursor *filename* option\n" if ((! defined $out_species_precursor) || (-e $out_species_precursor));
die "Need to specify a non existing non-target-species-mature output file via --outnonspeciesmature *filename* option\n" if ((! defined $out_nonspecies_mature) || (-e $out_nonspecies_mature));


my $kingdom	                        = &parse_organism($organism_file);
my ($precursor, undef) 			= &parse_fasta($precursor_file, $species, $kingdom);
my ($mature, $other_animal_mature)	= &parse_fasta($mature_file, $species, $kingdom);

my %assignment = (
    $out_species_mature    => $mature,
    $out_species_precursor => $precursor,
    $out_nonspecies_mature => $other_animal_mature,
    );

foreach my $outfile (keys %assignment)
{
    open(FH, ">", $outfile) || die "Unable to open file '$outfile': $!\n";
    foreach my $header (keys %{$assignment{$outfile}}){
	print FH "$header", "\n", $assignment{$outfile}{$header}, "\n";
    }
    close(FH) || die "Unable to close file '$outfile': $!\n";
}

sub parse_organism{
	my $po_file	= $_[0];
	my %po_metazoa;
	open(PO,"<",$po_file) || die;
	while(<PO>){
		chomp;
		next if(/^#/);
		my $po_line	= $_;
		my ($po_code,undef,undef,$po_kingdom,undef)	= split("\t",$po_line);
		if($po_kingdom =~ /^Metazoa/){
			$po_metazoa{$po_code} = "";
#			print "$po_code\n";
		}
	}
	close(PO) || die;
	return(\%po_metazoa);
}



sub parse_fasta{
	my $pf_file	= $_[0];
	my $pf_species	= $_[1];
	my %pf_kingdom	= %{$_[2]};
	my %pf_soi;	# species of interest
	my %pf_sni;	# species not of interest
	
	open(PF,"<",$pf_file) || die;
	my $pf_header;
	while(<PF>){
		chomp;
		my $pf_line	= $_;
		if(/^>/){
			$pf_header	= "";
			my @pf_split	= split("-",$pf_line);
			my $pf_code	= $pf_split[0];
			$pf_code	=~s/^>//;
			next unless (exists $pf_kingdom{$pf_code});
			$pf_header	= $pf_line;
		}
		else{
			next if($pf_header eq "");	
			my $pf_seq	= $pf_line;
			$pf_seq		=~ tr/Uu/Tt/;
			if($pf_header 	=~ /^>$species/i){
				$pf_soi{$pf_header}.=$pf_seq;
			}
			else{
				$pf_sni{$pf_header}.=$pf_seq;
			}
		}
	}
	close(PF) || die;
	return(\%pf_soi,\%pf_sni);
}



