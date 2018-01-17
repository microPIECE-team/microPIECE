#! /usr/bin/perl
use strict;
use warnings;
# read config file (fastq alignment files)

# output : csv len,rpm,type
# call 052_stacked_barplot.R

my $smrna_cfg	= $ARGV[0];
my %cfg_hash	= %{&read_cfg($smrna_cfg)};

foreach	my $condition	(keys %cfg_hash){
	my %replicates		= %{$cfg_hash{$condition}};
	my $count_replicates	= keys (%replicates);
	print "$condition : $count_replicates\n";
	foreach my $replicate	(keys %replicates){
		print "$replicate\n";
		my %alignment_fq	= %{$replicates{$replicate}};
		foreach my $fq_file (keys %alignment_fq){
			print "--> $fq_file : $alignment_fq{$fq_file}\n";
		}

	}

}






##########################################################################################
sub read_cfg{
	my $rc_file	= $_[0];
	my %rc_hash;	# {condition}{replicate}=path/to/file.fq
	open(RC,"<",$rc_file) || die;
	while(<RC>){
		chomp;
		next if(/^#/);
		next if(/^$/);
		my $rc_line	= $_;
		my ($rc_path,$rc_condition,$rc_replicate)	= split(" ",$rc_line);
#		my $rc_genome_aln	= 	$rc_path."_genome_aln.fq";
#		my $rc_ncRNA_aln	=	$rc_path."_ncRNA_aln.fq";
#		my $rc_miRNA_aln	=	$rc_path."_miRNA_aln.fq";
#		my $rc_mRNA_aln		= 	$rc_path."_mRNA_aln.fq";
#		my $rc_mRNA_ual		=	$rc_path."_mRNA_ual.fq";
		$rc_hash{$rc_condition}{$rc_path}{$rc_path."_genome_aln.fq"}	= 0;
                $rc_hash{$rc_condition}{$rc_path}{$rc_path."_ncRNA_aln.fq"}     = 0;
                $rc_hash{$rc_condition}{$rc_path}{$rc_path."_miRNA_aln.fq"}     = 0;
                $rc_hash{$rc_condition}{$rc_path}{$rc_path."_mRNA_aln.fq"}	= 0;
                $rc_hash{$rc_condition}{$rc_path}{$rc_path."_mRNA_ual.fq"}	= 0;

	}
	close(RC) || die;
	return(\%rc_hash);
}



sub count_fastq{
	

}
