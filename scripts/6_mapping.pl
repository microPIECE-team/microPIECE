#! /usr/bin/perl
use strict;
use warnings;

use File::Temp qw(tmpnam);


my $mir_file		= $ARGV[0];	#	tca_miRNA_mature_dna_novel.fa
my $target_seqs		= $ARGV[1];	#	map_clip_needle_out_max100_min"$i"_targetseqs.fa
my $prediction_out	= $ARGV[2];	#	map_clip_needle_out_max100_min"$i"_tarpredict_miranda.out

my $tmp_mir		= tmpnam();
my $tmp_miranda_out	= tmpnam();

#no thresholds defined so far

my %mir_hash	= %{&fasta_parser($mir_file)};



my %tar_hash;	# {tarID}	= \@(\@(mirnadaout),..)

foreach my $mir_header (keys %mir_hash){
	my $mir_seq 	= $mir_hash{$mir_header};
	open(MIR,">",$tmp_mir) || die;
	print MIR "$mir_header\n$mir_seq";
	close(MIR) || die;

	system("miranda $tmp_mir $target_seqs -quiet -out $tmp_miranda_out");
	open(IN,"<",$tmp_miranda_out) || die;
	while(<IN>){
		chomp;
		next unless (/^>[^>]/);
		my $in_line	= $_;
		my @in_array	= split(" ",$in_line);
		if(not exists $tar_hash{$in_array[1]}){
			my @tmp_array;
			push(@tmp_array,\@in_array);
			$tar_hash{$in_array[1]}	= \@tmp_array;
		}	
		else{
			my @tmp_array2	= @{$tar_hash{$in_array[1]}};
			push(@tmp_array2,\@in_array);
			$tar_hash{$in_array[1]} = \@tmp_array2;
		}
	}
	close(IN) || die;
}

#output
open (OUT,">",$prediction_out) || die;
print OUT "query_miRNA\ttarget_mRNA\tscore\tkcal/mol\tQuery-aln_start\tQuery-aln_end\tSubject-aln_start\tSubject_aln-End\tAl-Len\tSubject-Identity\tQuery-Identity\n";
foreach my $th_key (keys %tar_hash){
	my @th_array	= @{$tar_hash{$th_key}};
	foreach(@th_array){
		print OUT "@{$_}\n";
	}
}
close(OUT) || die;









#####################################################################################
sub fasta_parser{
        my $fp_file             = $_[0];
        my %fp_hash;
        my $fp_header           = "";
        open(FP,"<",$fp_file) || die;
        while(<FP>){
                chomp;
                my $fp_line = $_;
                if(/^>/){
                        my @fp_split    = split(" ",$fp_line);
                        $fp_header      = $fp_split[0];
                        if(not exists $fp_hash{$fp_header}){
                                $fp_hash{$fp_header} = "";
                        }
                        else{ 
                                print STDERR "$fp_header : has double entry\n";
                        }
                }
                else{
                        $fp_hash{$fp_header}    .= $fp_line;
                }
        }
        close(FP) || die;
        return(\%fp_hash);
}

