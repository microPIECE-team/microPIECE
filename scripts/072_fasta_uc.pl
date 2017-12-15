#! /usr/bin/perl
use strict;
use warnings;

my $fasta_file	= $ARGV[0];
my $counter	= 0;
open(FA,"<",$fasta_file) || die;
while(<FA>){
	chomp;
	my $fa_line	= $_;
	if(/^>/){
		print "$fa_line-$counter\n";
		$counter += 1;
	}
	else{
		my $fa_line	= $_;
		$fa_line	= uc($fa_line);
		print "$fa_line\n";
	}
}
close(FA) || die;
