#! /usr/bin/perl
use strict;
use warnings;


my $input_list		= $ARGV[0];	# comma separated replicates
my $condition_ID	= $ARGV[1];	# from db

my @input_array		= split(",",$input_list);
my %condition_hash;	#{condition_file}	= \%rep_hash;
foreach my $input_file	(@input_array){
	my %rep_hash;	 #{miRBaseID;mism;add;t5;t3;ambiguity} = freq/ambiguous
	open(TMP,"<",$input_file) || die;
	my $tmp_firstline	= <TMP>;
	my $tmp_total_read_count= 0;
	while(<TMP>){
		chomp;
		my @tmp_split	= split("\t",$_);
		my $tmp_seq	= $tmp_split[0];
		my $tmp_freq	= $tmp_split[2];
		my $tmp_id	= $tmp_split[3];
		my $tmp_mism	= $tmp_split[6];
		my $tmp_add	= $tmp_split[7];
		my $tmp_t5	= $tmp_split[8];
		my $tmp_t3	= $tmp_split[9];
		my $tmp_amb	= $tmp_split[14];
		my $tmp_div_c	= $tmp_freq/$tmp_amb;
	
		if($tmp_mism eq 0){
			$tmp_mism = 'NULL';
		}
		if($tmp_add eq 0){
			$tmp_add = 'NULL';
		}
		if($tmp_t5 eq 0){
			$tmp_t5 = 'NULL';
		}
		if($tmp_t3 eq 0){
			$tmp_t3 = 'NULL';
		}
		my $tmp_key	= "$tmp_id;$tmp_mism;$tmp_add;$tmp_t5;$tmp_t3;$tmp_seq";
		$tmp_total_read_count += ($tmp_freq/$tmp_amb);	# count all reads, divided by ambiguity, but also the non-isoforms
		next if(($tmp_mism eq 'NULL') && ($tmp_add eq 'NULL') && ($tmp_t5 eq 'NULL') && ($tmp_t3 eq 'NULL'));
		$rep_hash{$tmp_key}= $tmp_div_c;
	#	print ">>> $tmp_key\n";
	}
	close(TMP);
	# calc RPM
	my %rpm_hash;
	foreach my $rep_key (keys %rep_hash){
		my $rep_freq 	= $rep_hash{$rep_key};
		my $rep_rpm	= ($rep_freq/$tmp_total_read_count)*1000000;
		$rpm_hash{$rep_key}=$rep_rpm;
	}
	$condition_hash{$input_file}	= \%rpm_hash;

}
my %total_hash; #{miRBaseID;mism;add;t5;t3;ambiguity} = (RPM_con1 + RPM_con2 +...) / number of conditions 
my $condition_number = scalar(@input_array); # number of conditions : e.g. 4
foreach my $con_file (keys %condition_hash){
	my %tmp_rep_hash	= %{$condition_hash{$con_file}};
	foreach my $tmp_rep_key	(keys %tmp_rep_hash){
#		print"TMP REP KEY : $tmp_rep_key\n";
		if(not exists $total_hash{$tmp_rep_key}){
			$total_hash{$tmp_rep_key} = $tmp_rep_hash{$tmp_rep_key};
		}
		else{
			$total_hash{$tmp_rep_key} += $tmp_rep_hash{$tmp_rep_key};
		}
	}
}

# now divide every entry by the number of conditions
foreach my $total_key (keys %total_hash){
	my $total_hash_freq = $total_hash{$total_key}/$condition_number;
	my @total_key_split = split(";",$total_key);
#	print "@total_key_split\n";
#	print ">>>>>>>>>>>>>>>>> $total_key_split[0]\n";
	my $total_mir_ID = $total_key_split[0];
	print STDOUT "$total_mir_ID;$total_key_split[1];$total_key_split[2];$total_key_split[3];$total_key_split[4];$total_key_split[5];$total_hash_freq;$condition_ID\n";
}





