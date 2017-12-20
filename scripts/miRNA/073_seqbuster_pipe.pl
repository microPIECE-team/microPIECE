#! /usr/bin/perl
use strict;
use warnings;
my $fastq_path		= $ARGV[0];	# input folder
my $db			= $ARGV[1];	# database folder with hairpin.fa and miRNA.str
my $output		= $ARGV[2];	# output folder
my $species		= $ARGV[3];	# tca (like in mirbase)


my $collapse            = "/opt/anaconda2/bin/seqcluster collapse";
my $miraligner          = "java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s $species -freq";

opendir DIR, $fastq_path || die;
my @fastq_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR || die ;
$fastq_path	=~s/\/$//;

open(OUT,">","time_stamps.csv") || die;
foreach(@fastq_files){
	my $fastq_data	= $_;
	my $fastq_file	= "$fastq_path/$_";
	
	# collapse reads
	system("$collapse -f $fastq_file -o $output");
	
	# miraligner
	my $st = time();
	my $collapsed_file	= $fastq_data;
	$collapsed_file		=~s/\.fastq$//;
	$collapsed_file		=~s/\.fq$//;
	$collapsed_file		= "$output"."$collapsed_file"."_trimmed.fastq";
	my $aligned_file	= "$output"."$fastq_data";
	$aligned_file		=~s/\.fastq$//;
	$aligned_file		=~s/\.fq$//;
	system("$miraligner -i $collapsed_file -db $db -o $aligned_file");
	print OUT "TIME TAKEN $fastq_data:\t" . (time() - $st) . "\n";

}
close(OUT) || die;
