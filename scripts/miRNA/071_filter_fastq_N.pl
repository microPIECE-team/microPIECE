#! /usr/bin/perl
use strict;
use warnings;


my $fastq_file	= $ARGV[0];
my %fastq_hash;
my $fastq_header;
my $fastq_count	= 0;
open(FQ,"<",$fastq_file) || die;
while(<FQ>){
	chomp;
	my $fastq_line	 = $_;
	$fastq_count += 1;	
	if($fastq_count == 1){
		$fastq_header			= $fastq_line;
		my @fastq_tmp1			= ($fastq_line);
		$fastq_hash{$fastq_header} 	= \@fastq_tmp1;
	}	
	else{
		my @fastq_tmp2			= @{$fastq_hash{$fastq_header}};
		push(@fastq_tmp2,$fastq_line);
		$fastq_hash{$fastq_header}	= \@fastq_tmp2;
	}
	if($fastq_count == 4){
		$fastq_count = 0;
	}
}
close(FQ) || die;


foreach(keys %fastq_hash){
	my $fastq_key	= $_;
	my @fastq_array	= @{$fastq_hash{$fastq_key}};
	my $fastq_n	= "N";
	next if ( index ($fastq_array[1],$fastq_n) != -1);
	foreach(@fastq_array){
		print "$_\n";
	}
}
