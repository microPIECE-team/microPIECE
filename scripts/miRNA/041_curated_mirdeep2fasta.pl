#! /usr/bin/perl
use strict;
use warnings;
# updated score from , to . to avoid later processing bugs when transforming the sam file to a sam file with one line per hit
my $csv_file	= $ARGV[0];
my $class	= $ARGV[1]; # high_conf | low_conf

my $mature_file	= $csv_file;
$mature_file	=~ s/\.csv$/-mature.fa/;
my $hairpin_file= $csv_file;
$hairpin_file	=~s /\.csv$/-hairpin.fa/;



open(MATURE,">",$mature_file);
open(HAIRPIN,">",$hairpin_file);

open(CSV,"<",$csv_file);
while(<CSV>){
	chomp;
	my $csv_line	= $_;
	my @csv_split	= split(" ",$csv_line);
	my $source_tool	= $csv_split[0];
	my $tmp_id	= $csv_split[1];
	my $score	= $csv_split[2];
	$score		=~s/,/\./;
	my $ref_mir	= $csv_split[6];
	my $mature5p	= uc($csv_split[7]);
	my $mature3p	= uc($csv_split[8]);
	my $hairpin	= uc($csv_split[9]);
	my $pos		= $csv_split[10];
	
	$mature5p	=~ s/U/T/g;
	$mature3p	=~ s/U/T/g;
	$hairpin	=~ s/U/T/g;


	my $header	= ">tca-$pos|$score|$source_tool|$ref_mir|$class";

	print HAIRPIN "$header\n$hairpin\n";
	print MATURE "$header-5p\n$mature5p\n$header-3p\n$mature3p\n";
}
close(CSV);

close(MATURE);
close(HAIRPIN);
