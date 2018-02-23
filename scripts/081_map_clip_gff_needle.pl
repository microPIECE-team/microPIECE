#! /usr/bin/perl
use strict;
use warnings;

use File::Temp qw(tmpnam);

my $gff_csv	= $ARGV[0];	# use to transfer XM_AAE to XP_AAE
my $ortho_file	= $ARGV[1];	# use to transfer XP_AAE to XP_ortholog_in_TCA
my $clip_file	= $ARGV[2];	# use to check for clip ind XM_AAE
my $gff_csv2	= $ARGV[3];	# use to check for XM_TCA to XP_TCA
my $trans_file	= $ARGV[4];	# use to parse XM_TCA sequence for needle alignment
my $outfile     = $ARGV[5];

my $clip_len_cutoff=55;
# OUTPUT :
# XM_AAE	XP_AAE		#peaks		XM_TCA		XP_TCA		# peaks

my %gff_hash	= %{&GFF_parser($gff_csv)};	# {xp}			= xm
my %ortho_hash	= %{&ortho_parser($ortho_file)};# {XP_AAE1}		= \@(XP_TCA1,XP_TCA2,..)
my %clip_hash	= %{&fasta_parser($clip_file)};	# {header-count}	= seq
my %gff2_hash	= %{&GFF_parser($gff_csv2)};	# {XP_TCA1}		= XM_TCA1
my %trans_hash	= %{&fasta_parser($trans_file)};# {header}		= seq

open(OUT,">",$outfile) or die "could not write output\n";
print OUT "XM_QRY;XP_QRY\t||\tXM_TAR;XP_TAR\t||\tqry_header tar_header identity coverage QRY_len TAR_len qry_gaps tar_gaps total_gaps MM matches start stop seq\n";

my $tmp_mrna 	= tmpnam();
my $tmp_clip 	= tmpnam();
my $tmp_needle	= tmpnam();

foreach(sort keys %gff_hash){
	my $gff_xp_key		= $_;
	my $gff_xm_key		= $gff_hash{$gff_xp_key};
	next unless (exists $ortho_hash{$gff_xp_key});
	# check how many clip sequences map to mRNA XM	
	my @original_clip_seqs;	# from the original dataset, collect headers
	foreach(keys %clip_hash){
		my $clip_key	= $_;		# XM_...-counter
		my $clip_trim	= $_;
		my $clip_seq	= $clip_hash{$clip_key};
		my $clip_seq_len= length($clip_seq);
		next if ($clip_seq_len > $clip_len_cutoff);
		$clip_trim	=~s/-.*$//;
		if($clip_trim eq $gff_xm_key){
			push(@original_clip_seqs,$clip_key);
		}
	}
#	print "@original_clip_seqs\n";
	next unless (scalar(@original_clip_seqs) > 0);
	my @ortho_array	= @{$ortho_hash{$gff_xp_key}};	# XP_AME...
	foreach(@ortho_array){
		my $oa_xp_id	= $_;
		next unless (exists $gff2_hash{$oa_xp_id});
		my $oa_xm_id	= $gff2_hash{$oa_xp_id};
		open(MRNA,">",$tmp_mrna) or die;
		print MRNA ">$oa_xm_id\n$trans_hash{$oa_xm_id}\n";
		close(MRNA);
		foreach my $ocs_header (@original_clip_seqs){
#			print STDERR ">>>>: @original_clip_seqs\n";
#			my $ocs_header 	= $_;
			open(CLIP,">",$tmp_clip) or die;
			print CLIP ">$ocs_header\n$clip_hash{$ocs_header}";
			close(CLIP);
			system("needle -asequence $tmp_clip -bsequence $tmp_mrna  -datafile EDNACUSTOM -endweight Y -gapopen 5 -gapextend 2 -auto -aformat markx3 -outfile $tmp_needle");
			
			# parse needle output
			# (identity,coverage,gaps,mm, tar_start, tar_stop) 
			my $clip_seq_len	= length($clip_hash{$ocs_header});
			my @needle_array	= @{&parse_needle($tmp_needle,$clip_seq_len)};
			next unless (@needle_array);
			# ($pn_qry_header,$pn_tar_header,sprintf("%.3f",$pn_ident),sprintf("%.3f",$pn_coverage),$pn_qry_len,$pn_tar_seq_len,$pn_gaps,$pn_mm,$pn_match_nt,$pn_start,$pn_stop,$pn_tar_seq_clean)
			print OUT "$gff_xm_key;$gff_xp_key\t||\t$oa_xp_id;$oa_xm_id\t||\t@needle_array\n"; 
			
		}
		
	}
	# parse XP to XM for TCA as soon as csv of gff exists 
	#XM_original  ; XP_original	(clipheader-1,clipheader-2,..)	
}
close(OUT) || die;

################################################################################################
################################################################################################

sub parse_needle{
	my $pn_file	= $_[0];
	my $pn_qry_len	= $_[1];
	my $pn_ident	= 0;		# matching nts / qry_len * 100 
	my $pn_coverage	= 0;		
	my $pn_qry_gaps	= 0;		# count internal query gaps
	my $pn_tar_gaps	= 0;		# count internal target gaps
	my $pn_tot_gaps	= 0;		# total gaps = sum of qry and tar gaps
	my $pn_mm	= 0;		# count internal mm
	my $pn_start	= 0;		# start pos on target
	my $pn_stop	= 0;		# stop pos on target

	my $pn_match_nt		= 0;
	my $pn_aln_len;	

	my @pn_array;		#(pn_ident, pn_coverage, pn_gaps, pn_mm, pn_start, pn_stop, pn_tarseq)
	my $pn_count_header	= 0;
	my $pn_qry_header	= "";
	my $pn_qry_seq		= "";
	my $pn_tar_header	= "";
	my $pn_tar_seq		= "";

	open(PN,"<",$pn_file) || die;
	while(<PN>){
		chomp;
		my $pn_line	= $_;
		next if (/^#/);
		next if (/^$/);
		if (/^>/){
			if($pn_count_header == 0){
				$pn_qry_header = $pn_line;
				$pn_count_header += 1;
			}
			else{
				$pn_tar_header = $pn_line;
				$pn_count_header += 1;
			}
		}
		else{
			if($pn_count_header == 1){
				$pn_qry_seq .= $pn_line;
			}
			else{
				$pn_tar_seq .= $pn_line;
			}
		}
	}
	close(PN) || die;

		

	print "$pn_qry_header\t:\t$pn_qry_seq\n";
	print "$pn_tar_header\t:\t$pn_tar_seq\n";
	#################################################	
	#cut off 5' and 3' - 
	# start stop pos 
	my $pn_aln_len_seq	= $pn_qry_seq;
	$pn_aln_len_seq		=~ s/(^-*)//;
	$pn_start		= length($1);
	$pn_aln_len_seq		=~ s/(-*$)//;
	$pn_aln_len		= length($pn_aln_len_seq);	
	$pn_stop		= $pn_start+$pn_aln_len-1;
	#COVERAGE
#	my $pn_aln_len_seq_clean= $pn_aln_len_seq;
#	$pn_aln_len_seq_clean	=~s/-//g:
#	$pn_coverage	= $pn_aln_len_clean/$pn_qry_len*100;
	#################################################
	# QRY GAPS
	my @pn_aln_len_seq_split= split("",$pn_aln_len_seq);
	foreach(@pn_aln_len_seq_split){
		if($_ eq "-"){
			$pn_qry_gaps += 1;
		}
	}
	#################################################

	my @pn_qry_seq_split	= split("",$pn_aln_len_seq);
	#################################################
	# GET SUBSTRING OF TARGET SEQUENCE
	my @pn_tar_seq_split	= split("",$pn_tar_seq);
	return [] if (($pn_tar_seq_split[0] eq "-") or ($pn_tar_seq_split[-1] eq "-"));
	my $pn_tar_start	= -1;
	my $pn_tar_stop		= -1;
	my $pn_tar_pos		= 0;
	for(my $i = 0;$i<@pn_tar_seq_split;$i++){
		if($pn_tar_seq_split[$i] ne "-"){
			$pn_tar_pos ++;
		}
		if($i==$pn_start){
			$pn_tar_start = $pn_tar_pos;
		}
		if($i==$pn_stop){
			$pn_tar_stop = $pn_tar_pos;
		}	
	}


	my $pn_len		= $pn_stop - $pn_start + 1;
	my $pn_tar_seq_sub	= substr($pn_tar_seq, $pn_start, $pn_len);
	my $pn_tar_seq_len	= length($pn_tar_seq);
	my @pn_tar_seq_sub_split= split("",$pn_tar_seq_sub);
	my $pn_tar_seq_trim	= $pn_tar_seq_sub;
#	$pn_tar_seq_trim	=~ s/(^-*)//;
#	$pn_tar_seq_trim	=~ s/(-*$)//;
	my @pn_tar_seq_trim_split=split("",$pn_tar_seq_trim);
	# TAR GAPS
	foreach(@pn_tar_seq_trim_split){
		if($_ eq "-"){
			$pn_tar_gaps +=1;
		}
	}	
	# for output
	my $pn_tar_seq_clean	= $pn_tar_seq_sub;
	$pn_tar_seq_clean	=~s/-//g;
	#################################################
	# COMPARE MATCHES
	for(my $i=0;$i<@pn_qry_seq_split;$i++){
		my $pn_qry_char = $pn_qry_seq_split[$i];
		my $pn_tar_char = $pn_tar_seq_sub_split[$i];
		next if($pn_qry_char eq "-");
		next if($pn_tar_char eq "-");
		if($pn_qry_char ne $pn_tar_char){
			$pn_mm += 1;
		}
		else{
			$pn_match_nt	+= 1;
		}
	}

	$pn_stop	-= $pn_tar_gaps;	# correct stop position by internal gaps of target sequence
	$pn_ident	= $pn_match_nt 	/ ($pn_match_nt + $pn_mm) * 100;
	$pn_coverage	= ($pn_mm+$pn_match_nt)	/ $pn_qry_len * 100;
	$pn_tot_gaps	= $pn_tar_gaps + $pn_qry_gaps;
	@pn_array	= ($pn_qry_header,$pn_tar_header,sprintf("%.3f",$pn_ident),sprintf("%.3f",$pn_coverage),$pn_qry_len,$pn_tar_seq_len,$pn_qry_gaps,$pn_tar_gaps,$pn_tot_gaps,$pn_mm,$pn_match_nt,$pn_tar_start,$pn_tar_stop,$pn_tar_seq_clean);
#	print "@pn_array\n";
#	exit;

	return(\@pn_array);
}





# {XP_AAE1}	= \@(XP_TCA1,...)	
# {XP_right1}	= \@XP_left
# Species       Genes   Alg.-Conn.      GCF_000002335.3_Tcas5.2_protein.faa     GCF_000004015.4_AaegL3_protein.faa
sub ortho_parser{
	my $op_file	= $_[0];
	my %op_hash;
	open(OP,"<",$op_file) || die;
	while(<OP>){
		chomp;
		my $op_line		= $_;
		my @op_split		= split(" ",$op_line);
		my $op_left		= $op_split[3];
		my $op_right		= $op_split[4];
		my @op_left_split	= split(",",$op_left);
		my @op_right_split	= split(",",$op_right);
		foreach(@op_right_split){
			if(not exists $op_hash{$_}){
				$op_hash{$_}	= \@op_left_split;
			}
			else{
				print STDERR "$_ already exists!\n";
			}
		}
	}
	close(OP) || die;
	return(\%op_hash);
}

# {xp}	= xm
sub GFF_parser{
	my $gp_file	= $_[0];
	my %gp_hash;	# {xp}	= xm
	
	open(GP,"<",$gp_file) || die;
	while(<GP>){
		chomp;
		my $gp_line	= $_;
		next if(/^\t$/);
#		print "$gp_line\n";
		my @gp_split	= split(" ",$gp_line);
		my $gp_xm	= $gp_split[0];
		my $gp_xp	= $gp_split[1];
		$gp_hash{$gp_xp}= $gp_xm;
	}
	close(GP) || die;
	return(\%gp_hash);
}

# {id}  = fa_Seq
sub fasta_parser{
        my $fp_file             = $_[0];
        my %fp_hash;
        my $fp_header           = "";
	my @header              = ();
        my $fp_count            = 0;
        open(FP,"<",$fp_file) || die;
        while(<FP>){
                chomp;
                my $fp_line = $_;
                if(/^>/){
                        my @fp_split    = split(" ",$fp_line);
                        $fp_header      = $fp_split[0];
                        $fp_header      =~s/^>//;

			if(not exists $fp_hash{$fp_header}){
			    $fp_hash{$fp_header} = "";
			}
			else{	# fasta clip headers were unified --> should not be processed anymore
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

