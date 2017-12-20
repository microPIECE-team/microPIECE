#! /usr/bin/perl
# Author : Daniel Amsel - daniel.amsel@ime.fraunhofer.de
# This script performs a blast search and appends the query sequence length and the coverage at the end of the blast result
# Version 0.5
# modified output -> more values
# -------------------------------
# Version 0.4
# -p now specifies the blast
# type dynamically
# including help file
# -------------------------------
# Version 0.3
# query can now consist of more
# than one sequence line per 
# header
# -------------------------------
# Version 0.2
# added threads
# added print information
# -------------------------------
# Version 0.1
# first version
# -------------------------------
use strict;
use warnings;
use Getopt::Long;

if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
help();
exit;
}




my $blast;
my $query;
my $db;
my $out;
my $threads=1;


GetOptions      (       'p=s'			=>	\$blast			, 
			'query=s'       	=>      \$query 		,
                        'db=s'          	=>      \$db   			,
                        'out=s'         	=>      \$out   		,
                        'num_threads=i' 	=>      \$threads		);	


print "reading query file...";
my %qry_hash	= %{&fasta_to_hash($query)};
print "finished!\n";

print "reading db file...";
my %db_hash	= %{&fasta_to_hash($db)};
print "finished!\n";


my %qry_seq_hash 	= %{&fasta_to_seq($query)};

my %subject_seq_hash 	= %{&fasta_to_seq($db)};


print "Blasting...";
system("$blast -db $db -query $query -outfmt 6 -out $out -num_threads $threads -evalue 1000 -dust no -soft_masking false");
print "finished\n";
print "Adding sequence length and identity to output...";
open(OUT,">","$out.custom");
open(TMP,"<",$out);
while(<TMP>){
        chomp;
        my $line = $_;
        my @split_line  = split(' ',$line);     # ame-bantam      jcf7180018464735        100.00  20      0       0       1       20      54      35      0.003   40.1
        my $qry_header  = $split_line[0];
	my $db_header	= $split_line[1];
        my $hit_len     = $split_line[3];
	my $score	= $split_line[11];
	my $evalue	= $split_line[10];
        my $qry_seq_len = $qry_hash{$qry_header};
	my $db_seq_len  = $db_hash{$db_header};
	my $qry_coverage;
	my $db_coverage;
	if($blast eq "tblastn"){
		$qry_coverage	= $hit_len/$qry_seq_len*100;
		$db_coverage	= $hit_len/$db_seq_len*3*100;       
	}
	elsif($blast eq "blastx"){
		$qry_coverage 	= $hit_len/$qry_seq_len*3*100;
		$db_coverage	= $hit_len/$db_seq_len*100;
	}
	else{
		$qry_coverage	= $hit_len/$qry_seq_len*100;
		$db_coverage	= $hit_len/$db_seq_len*100;
	}
	my $query_seq	= $qry_seq_hash{$qry_header};
	my $subject_seq = $subject_seq_hash{$db_header};

	print OUT "$line\t$qry_seq_len\t$db_seq_len\t$qry_coverage\t$db_coverage\n";
}
close(TMP);
close(OUT);
print "finished\n";
print "Results at $out.custom\n";


sub help{ 
	print "-p blast-type\n-query input_query\n-db blast_DB\n-out output_file\n-num_threads #CPUs\n";
}



sub fasta_to_hash{
	my $fth_fasta_file	= $_[0];
	my $fth_header;
	my %fth_hash;
	open(FAF,"<",$fth_fasta_file);
	while(<FAF>){
	        chomp;
	        if (/^>/){
	                my @fth_split_array = split(' ',$_);
	                $fth_header = $fth_split_array[0];
	
	                $fth_header =~ s/>//g;
	                $fth_hash{$fth_header}=0;
	        }
	        else{
	                my $fth_len = length($_);
	                $fth_hash{$fth_header}+=$fth_len;
	        }
	}
	close(FAF);
	return(\%fth_hash);

}


sub fasta_to_seq{
	my $fts_fasta_file	= $_[0];
	my $fts_header;
	my %fts_hash;
	open(FAF,"<",$fts_fasta_file);
	while(<FAF>){
		chomp;
		my $line	= $_;
		if(/^>/){
			my @fts_split_array = split(' ',$line);
			$fts_header = $fts_split_array[0];
			$fts_header =~ s/>//;
			$fts_hash{$fts_header}="";
		}
		else{
			my $fts_seq = $line;
			$fts_hash{$fts_header}.=$fts_seq;
		}

	}
	close(FAF);
	return(\%fts_hash);
}
