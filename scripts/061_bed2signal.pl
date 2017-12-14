#!/usr/bin/perl
use strict;
use warnings;
# takes ff bedfile and extracts signal positions
# signal strength defined by user


my $bed_file		= $ARGV[0];
my $signal_strength	= $ARGV[1];

open(BED,"<",$bed_file) || die;
while(<BED>){
	chomp;
	if($.==1){
		print "$_\n";
		next;	
	}
	my $bed_line	= $_;	# chr	start	stop	info	.	strand
	my ($bed_chr, $bed_start, $bed_stop, $bed_info, undef, $bed_strand)	= split("\t",$bed_line);
	my @bed_nt_info		= split(";",$bed_info);	# length=XY;counts=pos_1sum/pos1_rep1/pos2_rep2/.../pos1_repN/,pos2_sum/pos2_rep1/...
	my $bed_nt_len		= $bed_nt_info[0];
	my $bed_nt_cnt		= $bed_nt_info[1];		# count=5/1/1/0/1/1/1,1/1/0/0/0/0/0,...
	$bed_nt_cnt		=~s/^counts=//;			#       5/1/1/0/1/1/1,1/1/0/0/0/0/0,...
	my @bed_nt_pos_list	= split(",",$bed_nt_cnt);	#	@(5/1/1/0/1/1/1 , 1/1/0/0/0/0/0 , ...)
	# create sum array
	my $bed_pos_cnt		= $bed_start;
	my %bed_pos_tmp_hash;	#{tmpID}=\@(pos1,pos2,pos3,...)
	my $bed_pos_tmp_ID	= 0;
	foreach my $bed_nt_pos (@bed_nt_pos_list){
		my @bed_nt_pos_rep_list	= split("/",$bed_nt_pos);	#@(5,1,1,0,1,1,1)
		my $bed_nt_pos_rep_sum	= $bed_nt_pos_rep_list[0];
		
		#check for signal strength cutoff
		if($bed_nt_pos_rep_sum >= $signal_strength){
			push(@{$bed_pos_tmp_hash{$bed_pos_tmp_ID}},"$bed_pos_cnt;$bed_nt_pos");
		}
		else{
			$bed_pos_tmp_ID++;
		}
		$bed_pos_cnt++;	
	}
	foreach(keys %bed_pos_tmp_hash){
		my $bed_pos_tmp_hash_id		= $_;
		my @bed_pos_tmp_hash_array	= @{$bed_pos_tmp_hash{$bed_pos_tmp_hash_id}};
		my @bed_pos_tmp_hash_positions	= ();	# 1000,1001,1002,...
		my @bed_pos_tmp_hash_info	= ();	# 5/1/1/0/1/1/1,1/1/0/0/0/0/0,...
		foreach(@bed_pos_tmp_hash_array){
			my ($bed_pos_tmp_hash_array_pos, $bed_pos_tmp_hash_array_nt) = split(";",$_);
			push(@bed_pos_tmp_hash_positions,$bed_pos_tmp_hash_array_pos);
			push(@bed_pos_tmp_hash_info,$bed_pos_tmp_hash_array_nt);
		}
		my $bed_pos_tmp_hash_info_str		= join(",",@bed_pos_tmp_hash_info);
		my @bed_pos_tmp_hash_positions_sort 	= sort { $a <=> $b } @bed_pos_tmp_hash_positions;
		my $bed_pos_tmp_hash_positions_start	= $bed_pos_tmp_hash_positions_sort[0];
		my $bed_pos_tmp_hash_positions_stop	= $bed_pos_tmp_hash_positions_sort[-1];
		print "$bed_chr\t$bed_pos_tmp_hash_positions_start\t$bed_pos_tmp_hash_positions_stop\t$bed_pos_tmp_hash_info_str\t.\t$bed_strand\n";
	}
}
close(BED) || die;


# output:
# chr	start	stop	sum1,sum2,sum3,...	.	strand
