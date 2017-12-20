#! /usr/bin/perl
use strict;
use warnings;

my $fastq2fasta		= 'fastq2fasta.pl';
my $remove_ws		= 'remove_white_space_in_id.pl';
my $collapse_reads	= 'collapse_reads_md.pl';		# collapse_reads.pl reads.fa mmu
my $bwa_index		= 'bwa index';
my $bwa_aln		= 'bwa aln';
my $bwa_samse		= 'bwa samse';
my $bwa_mult		= './xa2multi.pl';
my $bwa_sam_conv	= 'bwa_sam_converter.pl';
my $mirdeep		= 'miRDeep2.pl';

my $dir			= $ARGV[0];
my $out			= $ARGV[1];
my $ref_genome		= $ARGV[2];
my $mature_ref_mir      = $ARGV[3];
my $mature_other_mir    = $ARGV[4];
my $hairpin_ref_mir     = $ARGV[5];
my $cpu			= $ARGV[6];

if (not $dir =~/\/$/){
	$dir .= "/";
}
if (not $out =~/\/$/){
	$out .= "/";
}

my $ref_genome_no_ws	= $ref_genome;
$ref_genome_no_ws      .= "_noWhitespace.fa";
system("$remove_ws $ref_genome > $ref_genome_no_ws");

my $mature_ref_mir_no_ws	= $mature_ref_mir;
$mature_ref_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $mature_ref_mir > $mature_ref_mir_no_ws");
my $mature_other_mir_no_ws	= $mature_other_mir;
$mature_other_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $mature_other_mir > $mature_other_mir_no_ws");
my $hairpin_ref_mir_no_ws	= $hairpin_ref_mir;
$hairpin_ref_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $hairpin_ref_mir > $hairpin_ref_mir_no_ws");

opendir DIR, $dir;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR;

# bwa index                                                                                                                                                                                                        
system("$bwa_index $ref_genome_no_ws");


foreach(@files){
	my $file	= $_;
	my $inpath	= $dir.$file;
	my $outpath	= $out.$file;
# fastq2fasta.pl fastq file > outputfile
	my $fasta_file	= $outpath;
	$fasta_file	=~s/\./_/g;
	$fasta_file    .= ".fasta";
	print "FASTQ2FASTA\n";
	system("$fastq2fasta $inpath > $fasta_file"); 	
# remove_white_space_in_id.pl
	my $fasta_no_whitespace	= $fasta_file;
	$fasta_no_whitespace	=~s/\.fasta$/_noWhitespace.fasta/;
	print "NO WHITESPACE\n";
	system("$remove_ws $fasta_file > $fasta_no_whitespace");
# collapse_reads.pl reads.fa tca
	my $fasta_collapse	= $fasta_no_whitespace;
	$fasta_collapse		=~s/\.fasta$/_collapse.fasta/;
	print "COLLAPSE\n";
	system("$collapse_reads $fasta_no_whitespace tca > $fasta_collapse");
# bwa alignment
# bwa aln -f out.sai -l 8 -n 1 -o 0 -e 0 -k 1 -t 20 ref.fa reads.fq
	my $bwa_sai	= $fasta_collapse; 
	$bwa_sai	=~s/\.fasta$/_bwa.sai/;	
	print "BWA ALN\n";
	system("$bwa_aln -f $bwa_sai -n 1 -o 0 -e 0 -k 1 -t $cpu $ref_genome_no_ws $fasta_collapse");
# bwa samse (sai to sam)
# bwa samse -f out.sam ref.fa bwaout.sai reads.fq
	my $bwa_sam	= $bwa_sai;
	$bwa_sam	=~s/\.sai$/.sam/;
	print "BWA SAMSE\n";
	system("$bwa_samse -f $bwa_sam $ref_genome_no_ws $bwa_sai $fasta_collapse");
# bwa xa2multi
	my $bwa2multi	= $bwa_sam;
	$bwa2multi	=~s/\.sam$/_multi.sam/;
	print "BWA MULTI\n";
	system("$bwa_mult $bwa_sam > $bwa2multi");
# bwa_sam_converter.pl -i mapped.sam -o collapsed_read_output_file -a mapping_file_in_arf_format
	my $bwa_arf	= $bwa2multi;
	$bwa_arf	=~s/multi\.sam/multi_sam\.arf/;
	print "BWA SAM 2 ARF\n";
	system("$bwa_sam_conv -i $bwa2multi -o $fasta_collapse -a $bwa_arf");
#mirdeep
#miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs -t Mouse 2>report.log

	my $mirdeep2_log	= $out."run.log";
	print "MIRDEEP2\n";
	system("$mirdeep $fasta_collapse $ref_genome_no_ws $bwa_arf $mature_ref_mir_no_ws $mature_other_mir_no_ws $hairpin_ref_mir_no_ws -P");	
}



