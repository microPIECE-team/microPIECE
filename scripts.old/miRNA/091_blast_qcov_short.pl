#! /usr/bin/perl
# Author : Daniel Amsel - daniel.amsel@ime.fraunhofer.de
use strict;
use warnings;
use Getopt::Long;

if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
help();
exit;
}




my $query;
my $db;
my $out;
my $threads=1;


GetOptions      (       'query=s'       	=>      \$query 		,
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


# -strand plus
print "Blasting...";
system("blastn -db $db -query $query -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -out $out -num_threads $threads -word_size 4 -evalue 10000 -strand plus");
print "finished\n";
print "Adding sequence length and identity to output...";
open(OUT,">","$out.custom");
open(TMP,"<",$out);
while(<TMP>){
        chomp;
        my $line = $_;
	my ($qry_header,$db_header,$identity,$hit_len,$mm,$gapopen,$qry_start,$qry_end,$tar_start,$tar_end,$eval,$score,$qry_seq,$tar_seq) = split(' ',$line);
	next unless ($hit_len >= 10);
	next unless (($qry_start==1) && ($tar_start==1));
	
	my @qry_seq_array	= split("",$qry_seq);
	my @tar_seq_array	= split("",$tar_seq);

	my $bool_keep		= 1;	# 1 keep ; 0 discard
	for (my $i = 0; $i<10; $i++){
		if($qry_seq_array[$i] ne $tar_seq_array[$i]){
			$bool_keep = 0;
		}
	}
	my $gap_count		= 0;
	my $mm_count		= 0;
	if($hit_len > 10){
		for( my $i = 10; $i < scalar(@qry_seq_array) ; $i++){
			if(($qry_seq_array[$i] eq "-") or ($tar_seq_array[$i] eq "-")){	
				$gap_count++;
			}
			elsif(($qry_seq_array[$i] ne $tar_seq_array[$i])){
				$mm_count++;
			}
		}
	}	
	if(($mm_count > 1) or ($gap_count > 1)){
		$bool_keep = 0;
	}
	elsif(($mm_count == 1) and ($gap_count == 1)){
		$bool_keep = 0;
	}
	
	next unless ($bool_keep);

        my $qry_seq_len = $qry_hash{$qry_header};
	my $db_seq_len  = $db_hash{$db_header};
	my $qry_coverage= $hit_len/$qry_seq_len*100;
	my $db_coverage	= $hit_len/$db_seq_len*100;
	my $query_seq	= $qry_seq_hash{$qry_header};
	my $subject_seq = $subject_seq_hash{$db_header};

	print OUT "$line\t$qry_seq_len\t$db_seq_len\t$qry_coverage\t$db_coverage\n";
}
close(TMP);
close(OUT);
print "finished\n";
print "Results at $out.custom\n";


sub help{ 
	print "-query input_query\n-db blast_DB\n-out output_file\n-num_threads #CPUs\n";
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
