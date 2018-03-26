#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $mirdeep_csv;
my $mature_fasta;
my $precursor_copies = "";

GetOptions(
	"mirdeep_out=s"		=>\$mirdeep_csv,
	"mature_fasta=s"	=>\$mature_fasta,
	"precursor_copies=s"	=>\$precursor_copies) || die;

# split the list of given precursors that have genomic copies
my @precursor_list	= split(",",$precursor_copies);
my %precursor_hash;	#{precursor} = 1
foreach my $precursor (@precursor_list){
	$precursor =~s/mir/miR/;
	$precursor =~s/-\d$//;
	$precursor_hash{$precursor} += 1;
}

my %mature_hash		= %{&parse_fasta($mature_fasta)};
my %pair_hash		= %{&get_pairs(\%mature_hash)};

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
	#check for miRBase annotations that have only one mature sequence
	my $precursor_name 	= $mature_name;
	$precursor_name		=~s/-[35]p$//;
	next unless (exists $pair_hash{$precursor_name});
	next if( $pair_hash{$precursor_name} == 2);
	$mature_seq	= uc($mature_seq);
	$star_seq	= uc($star_seq);
	$hairpin_seq	= uc($hairpin_seq);
	my $p5_seq;
	my $p3_seq;
	# as miRDeep2 reports mature and star sequences instead of 5p and 3p, it is necessary to identify the arm of the mature and star sequence
	if(exists $mature_hash{$mature_name}){
		my $hairpin_mid	= (length($hairpin_seq))/2;
		my $mature_idx  = index($hairpin_seq,$mature_seq);
		# if the mature sequence index is located on the right side of the hairpin, its defined as 3p and the star sequence is defined as 5p
		if($mature_idx >= $hairpin_mid){
			$p5_seq = $star_seq;
			$p3_seq = $mature_seq;
		}
		# if not, then the mature sequence is defined as 5p and the star sequence as 3p
		else{
			$p5_seq = $mature_seq;
			$p3_seq = $star_seq;
		}
		print ">$precursor_name-5p\n$p5_seq\n>$precursor_name-3p\n$p3_seq\n";
		# remove already identified mature sequences from total mature sequences that cointain only one arm in the hairpin
		delete($mature_hash{$mature_name});
	}
}
close(CSV) || die;


foreach my $mature_ID (keys %mature_hash){
	my $precursor_ID	= $mature_ID;
	my $mature_arm		= substr($mature_ID,-2);
	$precursor_ID		=~s/-.p$//;

	# check if microRNA was provided by the user to have genomic copies 
	if(exists $precursor_hash{$precursor_ID}){
		# loop over all copies and copy mature sequences with new header : tca-mir-3811c-3p => tca-mir-3811c-1-3p and tca-mir-3811c-2-3p
		for(my $i=1;$i<=$precursor_hash{$precursor_ID};$i++){
			print ">$precursor_ID-$i-$mature_arm\n$mature_hash{$mature_ID}\n";
		}
	}
	#print all "normal" mature microRNAs
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


sub get_pairs{
	my %gp_mature	= %{$_[0]}; 	# {tca-miR-21-3p} ..
	my %gp_pair;			# {tca-miR-21} = count
	foreach my $gp_mature_key (keys %gp_mature){
		$gp_mature_key	=~s/-.p$//;
		$gp_pair{$gp_mature_key}+=1;
	}
	return(\%gp_pair);
}














