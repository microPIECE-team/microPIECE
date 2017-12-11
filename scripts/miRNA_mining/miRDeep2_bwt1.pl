#! /usr/bin/perl
use strict;
use warnings;

#my $makeblastdb		= 'makeblastdb';
#my $blastn		= 'blastn';

my $fastq2fasta		= 'fastq2fasta.pl';
my $remove_ws		= 'remove_white_space_in_id.pl';
my $collapse_reads	= 'collapse_reads_md.pl';		# collapse_reads.pl reads.fa mmu
my $mapper		= 'mapper.pl';
my $bowtie_build	= 'bowtie-build';


my $mirdeep			= 'miRDeep2.pl';

my $dir				= $ARGV[0];
my $out				= $ARGV[1];
my $ref_genome			= $ARGV[2];
my $mature_ref_mir_file      	= $ARGV[3];
my $mature_other_mir_file    	= $ARGV[4];
my $hairpin_ref_mir_file     	= $ARGV[5];
#my $nc_rna_ref_file		= $ARGV[6];
my $cpu				= $ARGV[6];

if (not $dir =~/\/$/){
	$dir .= "/";
}
if (not $out =~/\/$/){
	$out .= "/";
}

my $ref_genome_no_ws	= $ref_genome;
$ref_genome_no_ws      .= "_noWhitespace.fa";
system("$remove_ws $ref_genome > $ref_genome_no_ws");

my $mature_ref_mir		= &RNA2DNA($mature_ref_mir_file);
my $mature_ref_mir_no_ws	= $mature_ref_mir;
$mature_ref_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $mature_ref_mir > $mature_ref_mir_no_ws");
my $mature_other_mir		= &RNA2DNA($mature_other_mir_file);
my $mature_other_mir_no_ws	= $mature_other_mir;
$mature_other_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $mature_other_mir > $mature_other_mir_no_ws");
my $hairpin_ref_mir		= &RNA2DNA($hairpin_ref_mir_file);
my $hairpin_ref_mir_no_ws	= $hairpin_ref_mir;
$hairpin_ref_mir_no_ws	       .= "_noWhitespace.fa";
system("$remove_ws $hairpin_ref_mir > $hairpin_ref_mir_no_ws");

opendir DIR, $dir;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR;


my $bowtie_index	= "$ref_genome_no_ws-BWT_INDX";
system("$bowtie_build $ref_genome_no_ws $bowtie_index");



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

	# mapper (BOWTIE)
	my $bowtie_arf	= $fasta_collapse;
	$bowtie_arf	=~s/\.fasta$/.arf/;
	print "MAPPER\n";
	print "$mapper $fasta_collapse -c -q -n -l 17 -p $bowtie_index -t $bowtie_arf\n";
	system("$mapper $fasta_collapse -c -q -n -l 17 -p $bowtie_index -t $bowtie_arf");





#mirdeep
#miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs -t Mouse 2>report.log

	my $mirdeep2_log	= $out."run.log";
	print "MIRDEEP2\n";
	
	system("$mirdeep $fasta_collapse $ref_genome_no_ws $bowtie_arf $mature_ref_mir_no_ws $mature_other_mir_no_ws $hairpin_ref_mir_no_ws -P");	
}




# discard all 100% hits from fasta file
sub filter_blast{
	my $fb_blast_file	= $_[0];
	my $fb_fasta_file	= $_[1];
	
	my $fb_fasta_filtered	= $fb_fasta_file;
	$fb_fasta_filtered	=~s/\.fasta$/_filtered.fasta/;

	my %fb_blast_hash;
	my %fb_fasta_hash;


	my $fb_fasta_id;
	open(FA,"<",$fb_fasta_file);
	while(<FA>){
		chomp;
		my $fb_fasta_line	= $_;
		if(/^>/){
			$fb_fasta_id	= $fb_fasta_line;
			$fb_fasta_id	=~s/>//;
		}
		else{
			$fb_fasta_hash{$fb_fasta_id}=$fb_fasta_line;
		}
	}
	close(FA);
	
	open(BLAST,"<",$fb_blast_file);
	while(<BLAST>){
		chomp;
		my $fb_blast_line	= $_;
		my @fb_blast_array	= split(" ",$fb_blast_line);
		if($fb_blast_array[2] >= 100){	
			$fb_blast_hash{$fb_blast_array[0]}=$fb_blast_array[3];	# {query_coverage} = aln_length
		}
	}
	close(BLAST);


	foreach(keys %fb_fasta_hash){
		my $fb_fasta_key	= $_;
		if(exists $fb_blast_hash{$fb_fasta_key}){
			my $fb_fasta_seq_len	= length($fb_fasta_hash{$fb_fasta_key});
			if ($fb_fasta_seq_len == $fb_blast_hash{$fb_fasta_key}){
				delete $fb_fasta_hash{$fb_fasta_key};
				#print "QRY-LEN: $fb_fasta_seq_len | ALN-LEN: $fb_blast_hash{$fb_fasta_key} | $fb_fasta_key\n";
			}
		}
	
	}
	

	open(OUT,">",$fb_fasta_filtered);
	foreach(keys %fb_fasta_hash){
		print OUT ">$_\n";
		print OUT "$fb_fasta_hash{$_}\n";
	}
	close(OUT);
	return($fb_fasta_filtered);
}


sub RNA2DNA{
        my $r2d_rna_file        = $_[0];
        my $r2d_dna_file        = $r2d_rna_file."_dna.fa";
        open(OUT,">",$r2d_dna_file);
        open(R2D,"<",$r2d_rna_file);
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
        close(R2D);
        close(OUT);
        return($r2d_dna_file);
}












