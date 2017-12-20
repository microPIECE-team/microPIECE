#! /usr/bin/perl
use strict;
use warnings;


my $mirdeep_csv_bwa             = $ARGV[0];     # bwa
my $mirdeep_csv_bwt		= $ARGV[1];     # bwt
my $mirdeep_cutoff              = $ARGV[2];     # 4?

# MIRDEEP2 OUTPUT CSV
#[0]  provisional id    
#[1]  miRDeep2 score    
#[2]  estimated probability that the miRNA is a true positive   
#[3]  rfam alert        
#[4]  total read count  
#[5]  mature read count 
#[6]  loop read count   
#[7]  star read count   
#[8]  significant randfold p-value      
#[9]  mature miRBase miRNA      
#[10] example miRBase miRNA with the same seed  
#[11] UCSC browser      
#[12] NCBI blastn       
#[13] consensus mature sequence 
#[14] consensus star sequence   
#[15] consensus precursor sequence      
#[16] precursor coordinate



my %mirdeep_bwa_hash	= %{&parse_mirdeep($mirdeep_csv_bwa,$mirdeep_cutoff)};
my %mirdeep_bwt_hash	= %{&parse_mirdeep($mirdeep_csv_bwt,$mirdeep_cutoff)};




foreach(keys %mirdeep_bwa_hash){
	my $bwa_key	= $_;
	my @bwa_array	= @{$mirdeep_bwa_hash{$bwa_key}};	# \@(\@line1,\@line2,..)
	if(not exists $mirdeep_bwt_hash{$bwa_key}){
		$mirdeep_bwt_hash{$bwa_key}	= \@bwa_array;
	}
	else{
		my @bwt_array	= @{$mirdeep_bwt_hash{$bwa_key}};
		foreach(@bwa_array){
			push(@bwt_array,$_);
		}
		$mirdeep_bwt_hash{$bwa_key}	= \@bwt_array;
	}
	

}

foreach(keys %mirdeep_bwt_hash){
	my $bwt_key	= $_;
	my @bwt_list	= @{$mirdeep_bwt_hash{$bwt_key}};
	foreach(@bwt_list){
		my @bwt_entry	= @{$_};
		print "@bwt_entry\n";
	}

}



######################################

sub parse_mirdeep{
	my $pm_file	= $_[0];
	my $pm_cutoff	= $_[1];
	my %pm_hash;
	my $pm_bool	= 0;
	open(PM,"<",$pm_file);
	while(<PM>){
		chomp;
                my $pm_line    = $_;
                if (/^$/){
                        $pm_bool = 0;
                        next;
                }
                if(/^provisional/){
                        $pm_bool = 1;
                        next;
                }
		if($pm_bool == 1){
			my @pm_split		= split("\t",$pm_line);
			my $pm_prov_id		= $pm_split[0];
			my $pm_score		= $pm_split[1];
			my $pm_mature_count	= $pm_split[5];
			my $pm_loop_count	= $pm_split[6];
			my $pm_star_count	= $pm_split[7];
			my $pm_mirbase_ref	= $pm_split[10];
			my $pm_mature_seq	= $pm_split[13];
			my $pm_star_seq		= $pm_split[14];
			my $pm_precursor_seq	= $pm_split[15];
			my $pm_precursor_coord	= $pm_split[16];	
#			next if($pm_score > 6);
#			next if($pm_score < 4);
			next if ($pm_score < $pm_cutoff);
#			next if ($pm_precursor_coord =~ /^NW/);
			my @pm_list	= ($pm_file,$pm_prov_id,$pm_score,$pm_mature_count,$pm_loop_count,$pm_star_count,$pm_mirbase_ref,$pm_mature_seq,$pm_star_seq,$pm_precursor_seq,$pm_precursor_coord);
			if(not exists $pm_hash{$pm_mature_seq}){
				my @pm_array	= ();
				push(@pm_array,\@pm_list);
				$pm_hash{$pm_mature_seq}	= \@pm_array;
			}
			else{
				my @pm_array2	= @{$pm_hash{$pm_mature_seq}};
				push(@pm_array2,\@pm_list);
				$pm_hash{$pm_mature_seq}	= \@pm_array2;
			}
		}
	}
	close(PM);	
	return(\%pm_hash);
}
