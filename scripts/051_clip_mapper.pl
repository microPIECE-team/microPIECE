#! /usr/bin/perl
use strict;
use warnings;

my $bed_file	= $ARGV[0];
my $gff_file	= $ARGV[1];

my %gff_hash	= %{&gff_parser($gff_file)};

my $genomic_loci	= 0;	# clip peaks inside genome regions
my $intergenic_loci	= 0;	# clip peaks in intergenic regions
my $putative_utr	= 0;	# clip peaks in putative UTR
my $peak_count		= 0;	# count peaks
open (BED,"<",$bed_file) || die;
while(<BED>){
	chomp;
	# ignore comment lines
	if ($_ =~ /^#/)
	{
	    print $_, "\n";
	    next;
	}
	my $bed_line	= $_;
	my @bed_split	= split("\t",$bed_line);
	my $bed_chr	= $bed_split[0];
	my $bed_start	= $bed_split[1];
	my $bed_stop	= $bed_split[2];
        my $bed_info    = $bed_split[3];
        my $bed_strand  = $bed_split[5];	
	my $bed_len	= $bed_stop-$bed_start+1; 	

	my $bed_center = int($bed_split[1]+$bed_len/2); # center of the bed peak

	# gff_hash contains array refs: element [1] = feature start, [2] feature stop, [4] strand
	# features which cover at least half of the bed region are extracted and joined by ,
	my $target_annotation = join(",", map { $_->[3] } grep { (($bed_center <= $_->[2]) && ($bed_center >= $_->[1])&&($bed_strand eq $_->[4])) } (@{$gff_hash{$bed_chr}}));

	print join("\t", ($bed_chr, $bed_start, $bed_stop, join(";", ($bed_info, "annotation=".$target_annotation)), ".", $bed_strand)), "\n";
}
close(BED) || die;


sub gff_parser{
	my $gp_file	= $_[0];
	my %gp_hash;	# {count}	= \@(chr,start,stop,id)
	my $gp_count	= 0;
	open(GP,"<",$gp_file)|| die ;
	while(<GP>){
		next if (/^#/);
		chomp;
		my $gp_line	= $_;
		my @gp_split	= split(" ",$gp_line);
		next unless ($gp_split[2] eq "mRNA");
		my $gp_info	= $gp_split[8];
		if($gp_info =~ /;Name=/){
			my @gp_info_split	= split(";Name=",$gp_info);
			my $gp_info_name	= $gp_info_split[1];
			$gp_info_name		=~ s/;.*//; #XM_...
#			print "$gp_info_name\n";	
			# CHR, start, stop , XM, strand
			my @gp_array	= ($gp_split[0],$gp_split[3],$gp_split[4],$gp_info_name,$gp_split[6]);
			push(@{$gp_hash{$gp_split[0]}},	\@gp_array);	# {chr} = \@( \@(mrna), \@(mrna2), ...)
#			$gp_count      += 1;
		}
		else{
			print STDERR  "NO NAME FOUND IN : $gp_info\n";
		}
	}
	close(GP) || die;
	return(\%gp_hash);
}
