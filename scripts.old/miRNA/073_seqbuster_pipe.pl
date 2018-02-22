#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $fastq_path;
my $db;
my $output;
my $species;

GetOptions(
	"fq_path=s"		=>\$fastq_path,
	"db_path=s"		=>\$db,
	"out_path=s"		=>\$output,
	"species=s"		=>\$species) || die;

my $collapse            = "/opt/anaconda2/bin/seqcluster collapse";
my $miraligner          = "java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s $species -freq";

opendir DIR, $fastq_path || die;
my @fastq_files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR || die ;
$fastq_path	=~s/\/$//;

foreach(@fastq_files){
	my $fastq_data	= $_;
	my $fastq_file	= "$fastq_path/$_";
	
	# collapse reads
	system("$collapse -f $fastq_file -o $output");
        if($? != 0){
                print STDERR "Execution fail : $!\n";
                die;
        }
	
	# miraligner
	my $collapsed_file	= $fastq_data;
	$collapsed_file		=~s/\.fastq$//;
	$collapsed_file		=~s/\.fq$//;
	$collapsed_file		= "$output"."$collapsed_file"."_trimmed.fastq";
	my $aligned_file	= "$output"."$fastq_data";
	$aligned_file		=~s/\.fastq$//;
	$aligned_file		=~s/\.fq$//;
	system("$miraligner -i $collapsed_file -db $db -o $aligned_file");
        if($? != 0){
                print STDERR "Execution fail : $!\n";
                die;
        }

}
