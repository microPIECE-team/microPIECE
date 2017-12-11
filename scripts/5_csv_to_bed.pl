#! /usr/bin/perl
use strict;
use warnings;

use File::Temp qw(tmpnam);


my $csv_file		= $ARGV[0];	#	map_clip_needle_out_max100_min"$i".csv
my $target_seqs		= $ARGV[1];	#	map_clip_needle_out_max100_min"$i"_targetseqs.fa


open(OUT,">",$target_seqs) || die;
open(CSV,"<",$csv_file) || die;

my @csv_out;

while(<CSV>){
	chomp;
	next if ($. == 1);
	my $csv_line	= $_;
	my @csv_array	= split(" ",$csv_line);	
	my $csv_qry	= $csv_array[0];
	my $csv_tar	= $csv_array[2];
	my $csv_clip	= $csv_array[4];
	my $csv_identity= $csv_array[8];
	my $csv_coverage= $csv_array[9];
	my $csv_qry_len	= $csv_array[10];
	my $csv_tar_len	= $csv_array[11];
	my $csv_qry_gaps= $csv_array[12];
	my $csv_tar_gaps= $csv_array[13];
	my $csv_tot_gaps= $csv_array[14];
	my $csv_MM	= $csv_array[15];
	my $csv_matches	= $csv_array[16];
	my $csv_start	= $csv_array[17];
	my $csv_stop	= $csv_array[18];
	my $csv_seq	= $csv_array[19];
	my $csv_header	= "QRY:$csv_qry|TAR:$csv_tar|CLIP:$csv_clip|IDENTITY:$csv_identity|COVERAGE:$csv_coverage|QRY_GAPS:$csv_qry_gaps|TAR_GAPS:$csv_tar_gaps|TOTAL_GAPS:$csv_tot_gaps|MM:$csv_MM|MATCHES:$csv_matches|TAR_START:$csv_start|TAR_STOP:$csv_stop";
	#print OUT "$csv_header\n$csv_seq\n";
	my (undef,$csv_tar_mRNA_ID) = split(";",$csv_tar);
	push(@csv_out,[$csv_tar_mRNA_ID,$csv_start,$csv_stop,$csv_header]);
	
}

@csv_out = sort{$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} (@csv_out);

foreach(@csv_out){
	print OUT join("\t",@$_),"\n";
}

close(CSV) || die;
close(OUT) || die;


