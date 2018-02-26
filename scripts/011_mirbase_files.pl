#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $precursor_file;
my $mature_file;
my $species;		# tca
my $organism_file;	# organisms.txt from mirbase
my $out;

GetOptions(
	"organisms=s"		=> \$organism_file,
	"species=s"		=> \$species,
        "precursor_file=s"      => \$precursor_file,
        "mature_file=s" 	=> \$mature_file,
	"out=s"			=> \$out) || die;




my %kingdom	= %{&parse_organism($organism_file)};


my ($tca_precursor_ref, undef) 			= &parse_fasta($precursor_file,$species,\%kingdom);
my ($tca_mature_ref,$other_animal_mature)	= &parse_fasta($mature_file,$species,\%kingdom);

my %tca_precursor				= %{$tca_precursor_ref};
my %tca_mature					= %{$tca_mature_ref};

my %other_animal_mature				= %{$other_animal_mature};

open(TMP,">","$out"."mature_mirbase.fa") || die;
foreach(keys %tca_mature){
	print TMP "$_\n$tca_mature{$_}\n";
}
close(TMP) || die;

open(TMP,">","$out"."precursor_mirbase.fa") || die;
foreach(keys %tca_precursor){
	print TMP "$_\n$tca_precursor{$_}\n";
}
close(TMP) || die;

open(TMP,">","$out"."mature.fa-no-speciesB.fa") || die;
foreach(keys %other_animal_mature){
	print TMP "$_\n$other_animal_mature{$_}\n";
}
close(TMP) || die;

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
			$pf_seq		=~ tr/U/T/;
			if($pf_header 	=~ /^>$species/){
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



