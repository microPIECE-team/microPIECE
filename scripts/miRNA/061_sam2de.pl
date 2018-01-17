#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# sam file to DE.csv file for R script 
my $cfg_file;
my $mature_mir_file;
GetOptions(
	"cfg=s"		=> \$cfg_file,
	"mature_file=s"	=> \$mature_mir_file) || die;

# output : rpm; condition; microRNA


my %sam_hash		= %{&read_cfg($cfg_file)};	#{rep_name}	= \@(sam1,sam2,sam3,sam4)
my %mir_hash		= %{&read_fasta($mature_mir_file)};	# {mirID} = 0;
my @mir_hash_keys	= keys %mir_hash;

print "rpm;condition;miRNA\n";

foreach(keys %sam_hash){	# {condition} = \@(rep1,rep2,rep3,rep4)
	my $sh_condition	= $_;
	my @sh_files		= @{$sam_hash{$sh_condition}};
	my $sh_rep_count	= scalar(@sh_files);		# number of replicates --> total mir count of all 4 rep / rep_count
	my %con_rpm_mir_hash	= %mir_hash;			# {mirID} = sum of RPM of all 4 replicates, divided by 4
	foreach(@sh_files){	# loop through replicates
		# parse sam	#{read}	= \@(target1,target2,..)
		my $sh_sam_file		= $_;				# replicate file 
		my %sh_read_hash	= %{&read_sam($sh_sam_file)};	# {read} = \@(tar1,tar2,tar3,..)
		my %sh_mir_hash		= %mir_hash;			# {mirID} =  RPM of each miR  of one replicate
		my $sh_read_total_count	= 0;
		foreach(keys %sh_read_hash){	# loop sam file of replicate
			my $sh_read_key		= $_;			# readID
			my @sh_read_tar_list	= @{$sh_read_hash{$sh_read_key}};	# \@(tar1,tar2,tar3,tar4)
		#	print ">>> @sh_read_tar_list\n";
			my $sh_tarcount		= scalar(@sh_read_tar_list);		# 4
			my $sh_multi_count	= 1/$sh_tarcount;		# 1 read 4 targets = 0.25 count per target
			foreach(@sh_read_tar_list){				# add the count to each miR ID 
				my $sh_read_tar_id	= $_;
				$sh_mir_hash{$sh_read_tar_id}	+= $sh_multi_count;	# per replicate
			}
			$sh_read_total_count	+= 1;			# count each read 
		}
		# loop through replicate mir hash and compute RPM of replicate then add it to condition_rpm_mir_hash
		foreach(keys %sh_mir_hash){	# per replicate
			my $sh_rep_mir_id	= $_;
			my $sh_rep_mir_count	= $sh_mir_hash{$sh_rep_mir_id};
			my $sh_rep_mir_RPM	= $sh_rep_mir_count / $sh_read_total_count * 1000000;
			$con_rpm_mir_hash{$sh_rep_mir_id}	+= $sh_rep_mir_RPM;
		#	print "$sh_rep_mir_id | mir count : $sh_rep_mir_count | mir PRM $sh_rep_mir_RPM\n";
		}	
	}
	foreach(keys %con_rpm_mir_hash){
		my $crm_mirID			= $_;
		my $crm_RPM_sum			= $con_rpm_mir_hash{$crm_mirID};
		my $crm_RPM_mean		= $crm_RPM_sum / $sh_rep_count;
		$con_rpm_mir_hash{$crm_mirID}	= $crm_RPM_mean;
	#	print "$crm_mirID | $crm_RPM_sum | $crm_RPM_mean\n";
	}
	foreach(@mir_hash_keys){
		print "$con_rpm_mir_hash{$_};$sh_condition;$_\n";
	}	
}

# {mirID}	= 0;
sub read_fasta{
	my $rf_file	= $_[0];
	my %rf_hash;	#{mirID} = 0
	open(RF,"<",$rf_file) || die;
	while(<RF>){
		chomp;
		next unless (/^>/);
		my $rf_line		= $_;
		my @rf_split		= split(" ",$rf_line);
		my $rf_mirID		= $rf_split[0];
		$rf_mirID		=~s/^>//;		# remove > 
		$rf_hash{$rf_mirID}	= 0;
	}
	close(RF) || die;
	return(\%rf_hash);
}





#	{read} = \@(tar1_ID,tar2_ID,..)
sub read_sam{
	my $rs_file	= $_[0];
	my %rs_hash;
	open(RS,"<",$rs_file) || die;
	while(<RS>){
		chomp;
		next if (/^@/);
		my $rs_line	= $_;
		my @rs_split	= split(" ",$rs_line);
		my $rs_readID	= $rs_split[0];
		my $rs_mirID	= $rs_split[2];
		if(not exists $rs_hash{$rs_readID}){
			my @rs_list1		= ();
			push(@rs_list1,$rs_mirID);
			$rs_hash{$rs_readID}	= \@rs_list1;
		}
		else{
			my @rs_list2		= @{$rs_hash{$rs_readID}};
			push(@rs_list2,$rs_mirID);
			$rs_hash{$rs_readID}	= \@rs_list2;
		}
	}
	close(RS) || die;
	return(\%rs_hash);
}




# {rep_name}	= \@(sam1,sam2,sam3,sam4)
sub read_cfg{
	my $rc_cfg_file	= $_[0];	# config file
	my %rc_hash;			# {rep_name}	= \@(sam1,sam2,sam3,sam4);
	open(RC,"<",$rc_cfg_file) || die;
	while(<RC>){
		chomp;
		next if(/^#/);
		my $rc_line		= $_;
		my @rc_split		= split(" ",$rc_line);
		my $rc_file		= $rc_split[0];
		my $rc_condition	= $rc_split[1];
#		push(@{$rc_hash{$rc_condition}},$rc_file);

		if(not exists $rc_hash{$rc_condition}){
			my @rc_list1		= ();
			push(@rc_list1,$rc_file);
			$rc_hash{$rc_condition}	= \@rc_list1;
		}
		else{
			my @rc_list2		= @{$rc_hash{$rc_condition}};
			push(@rc_list2,$rc_file);
			$rc_hash{$rc_condition}	= \@rc_list2;
		}

	}
	close(RC) || die;
	return(\%rc_hash);
}

