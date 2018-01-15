#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my $fastq2fasta		= 'fastq2fasta.pl';
my $remove_ws		= 'remove_white_space_in_id.pl';
my $collapse_reads	= 'collapse_reads_md.pl';		# collapse_reads.pl reads.fa mmu
my $mapper		= 'mapper.pl';
my $bowtie_build	= 'bowtie-build';


my $mirdeep			= 'miRDeep2.pl';

my $dir;
my $out;
my $ref_genome;
my $mature_ref_mir_file;
my $mature_other_mir_file;
my $hairpin_ref_mir_file;
my $cpu;

GetOptions(
    	"dir=s"  			=> \$dir,
    	"out=s"     			=> \$out,
    	"ref_genome=s" 			=> \$ref_genome,
    	"species_mature_miRs=s"		=> \$mature_ref_mir_file,
	"other_mature_miRs=s"		=> \$mature_other_mir_file,
	"species_precursor_mirs=s"	=> \$hairpin_ref_mir_file,
	"threads=i"			=> \$cpu) || die;


if (not $dir =~/\/$/){
	$dir .= "/";
}
if (not $out =~/\/$/){
	$out .= "/";
}

my $ref_genome_no_ws	= $ref_genome;
$ref_genome_no_ws      .= "_noWhitespace.fa";
system("$remove_ws $ref_genome > $ref_genome_no_ws");
if($? == -1){
	print STDERR "Execution fail : $!\n";
        die;
}



my $mature_ref_mir		= &RNA2DNA($mature_ref_mir_file);
my $mature_ref_mir_no_ws	= $mature_ref_mir;
$mature_ref_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $mature_ref_mir > $mature_ref_mir_no_ws");
if($? == -1){
        print STDERR "Execution fail : $!\n";
        die;
}

my $mature_other_mir		= &RNA2DNA($mature_other_mir_file);
my $mature_other_mir_no_ws	= $mature_other_mir;
$mature_other_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $mature_other_mir > $mature_other_mir_no_ws");
if($? == -1){
        print STDERR "Execution fail : $!\n";
        die;
}
my $hairpin_ref_mir		= &RNA2DNA($hairpin_ref_mir_file);
my $hairpin_ref_mir_no_ws	= $hairpin_ref_mir;
$hairpin_ref_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $hairpin_ref_mir > $hairpin_ref_mir_no_ws");
if($? == -1){
        print STDERR "Execution fail : $!\n";
        die;
}

opendir DIR, $dir;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR;


my $bowtie_index	= "$ref_genome_no_ws-BWT_INDX";
system("$bowtie_build $ref_genome_no_ws $bowtie_index");
if($? == -1){
        print STDERR "Execution fail : $!\n";
        die;
}



foreach(@files){
	my $file	= $_;
	my $inpath	= $dir.$file;
	my $outpath	= $out.$file;
	my $fasta_file	= $outpath;
	$fasta_file	=~s/\./_/g;
	$fasta_file    .= ".fasta";
	print "FASTQ2FASTA\n";
	system("$fastq2fasta $inpath > $fasta_file"); 
	if($? == -1){
	        print STDERR "Execution fail : $!\n";
       		die;
	}	
	my $fasta_no_whitespace	= $fasta_file;
	$fasta_no_whitespace	=~s/\.fasta$/_noWhitespace.fasta/;
	print "NO WHITESPACE\n";
	system("$remove_ws $fasta_file > $fasta_no_whitespace");
	if($? == -1){
	        print STDERR "Execution fail : $!\n";
	        die;
	}
	my $fasta_collapse	= $fasta_no_whitespace;
	$fasta_collapse		=~s/\.fasta$/_collapse.fasta/;
	print "COLLAPSE\n";
	system("$collapse_reads $fasta_no_whitespace tca > $fasta_collapse");
	if($? == -1){
       		print STDERR "Execution fail : $!\n";
	        die;
	}
	# mapper (BOWTIE)
	my $bowtie_arf	= $fasta_collapse;
	$bowtie_arf	=~s/\.fasta$/.arf/;
	print "MAPPER\n";
	system("$mapper $fasta_collapse -c -q -n -l 17 -p $bowtie_index -t $bowtie_arf");
	        if($? == -1){
                print STDERR "Execution fail : $!\n";
                die;
        }

	my $mirdeep2_log	= $out."run.log";
	print "MIRDEEP2\n";
	
	system("$mirdeep $fasta_collapse $ref_genome_no_ws $bowtie_arf $mature_ref_mir_no_ws $mature_other_mir_no_ws $hairpin_ref_mir_no_ws -P");
        if($? == -1){
                print STDERR "Execution fail : $!\n";
                die;
        }
	
}




sub RNA2DNA{
        my $r2d_rna_file        = $_[0];
        my $r2d_dna_file        = $r2d_rna_file."_dna.fa";
        open(OUT,">",$r2d_dna_file) || die;
        open(R2D,"<",$r2d_rna_file) || die;
        while(<R2D>){
                chomp;
                my $r2d_line    = $_;
                if(/^>/){
                        print OUT "$r2d_line\n";
                }
                else{
                        my $r2d_rna_seq = $_;
                        my $r2d_dna_seq = $r2d_rna_seq;
                        $r2d_dna_seq    =~s/U/T/g;
                        print OUT "$r2d_dna_seq\n";
                }
        }
        close(R2D) || die;
        close(OUT) || die;
        return($r2d_dna_file);
}












