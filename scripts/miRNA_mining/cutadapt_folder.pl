#! /usr/bin/perl
use strict;
use warnings;

my $cutadapt	= "cutadapt -a TGGAATTCTCGGGTGCCAAGG -g GTTCAGAGTTCTACAGTCCGACGATC --trim-n --minimum-length 17 --maximum-length 40"; # $path > $trimmed
my $infolder	= $ARGV[0];
my $outfolder	= $ARGV[1];

opendir DIR, $infolder;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR;


foreach(@files){
	my $file	= $_;
	my $inpath	= $infolder.$file;
	my $outpath	= $outfolder.$file;
	$outpath	=~ s/\.fastq$/_trim_17_40_N_5_3_adapter.fastq/;
	
	system("$cutadapt $inpath -o $outpath");

}


