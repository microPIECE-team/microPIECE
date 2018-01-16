#! /usr/bin/perl
use strict;
use warnings;
use GetOpt::Long;


GetOptions(
	"mirdeep_out=s"	=>\$mirdeep_csv,
	"mature_fasta=s"=>\$mature_fasta,
	"precursor_copies=s"=\$precursor_copies) || die;

my @precursor_list	= split(",",$precursor_copies);
my %precursor_hash;	#{precursor} = 1
foreach my $precursor (@precursor_list){
	$precursor =~s/mir/miR/;
	$precursor =~s/-\d$//;
	$precursor_hash{$precursor} += 1;
}

my %mature_hash		= %{&parse_fasta($mature_fasta)};


open(CSV,"<",$mirdeep_csv) || die;
my $line_bool = 0;
while(<CSV>){
	chomp;
	if(/^mature miRBase miRNAs detected by miRDeep2/){
		$line_bool = 1;
		next;
	}
	next if ($line_bool ==0);
	next if (/^tag/);
	if(/^$/){
		$line_bool = 0;
		next;
	}
	my (undef, undef, undef, undef, undef, undef, undef, undef, undef, $mature_name, undef, undef, undef, $mature_seq, $star_seq, $hairpin_seq, undef) = split("\t",$_);
	next if ($mature_name =~/.+p$/);
	$mature_seq	= uc($mature_seq);
	$star_seq	= uc($star_seq);
	$hairpin_seq	= uc($hairpin_seq);
	my $p5_seq;
	my $p3_seq;
	if(exists $mature_hash{$mature_name}){
		my $hairpin_mid	= (length($hairpin_seq))/2;
		my $mature_idx  = index($hairpin_seq,$mature_seq);
		if($mature_idx >= $hairpin_mid){
			$p5_seq = $star_seq;
			$p3_seq = $mature_seq;
		}
		else{
			$p5_seq = $mature_seq;
			$p3_seq = $star_seq;
		}
		print ">$mature_name-5p\n$p5_seq\n>$mature_name-3p\n$p3_seq\n";
		delete($mature_hash{$mature_name});
	}
}
close(CSV) || die;

foreach my $mature_ID (keys %mature_hash){
	my $precursor_ID	= $mature_ID;
	my $mature_arm		= substr($mature_ID,-2);
	$precursor_ID		=~s/-.p$//;
	
	if($precursor_hash{$precursor_ID}){
		for(my $i=1;$i<=$precursor_hash{$precursor_ID};$i++){
			print ">$precursor_ID-$i-$mature_arm\n$mature_hash{$mature_ID}\n";
		}
	}
	else{
		print ">$mature_ID\n$mature_hash{$mature_ID}\n";
	}
}




sub parse_fasta{
	my $pf_file	= $_[0];
	my %pf_hash;
	my $pf_header;
	open(PF,"<",$pf_file) || die;
	while(<PF>){
		chomp;
		if(/^>/){
			$pf_header = (split(" ",$_))[0];
			$pf_header =~s/^>//;
		}
		else{
			$pf_hash{$pf_header} = $_
		}
	}
	close(PF)|| die;
	return(\%pf_hash);
}
