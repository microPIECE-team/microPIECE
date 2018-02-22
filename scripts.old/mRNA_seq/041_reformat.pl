#! /usr/bin/perl
use strict;
use warnings;

# input:
# isoforms.fpkm_tracking from cuffdiff
# gff file

# output: rpm;condition;NCBI_mRNA_ID

my $cuff_file	= $ARGV[0];
my $gff_file	= $ARGV[1];


my %gff_hash; # {rnaID} = XM_..

open(GFF,"<",$gff_file) || die;
while(<GFF>){
        chomp;
        next if(/^#/);
        my $gff_line            = $_;
        my @gff_split           = split("\t",$gff_line);
        next unless (($gff_split[2] eq "mRNA") or ($gff_split[2] eq "transcript"));
        my $gff_info            = $gff_split[8];
        my @gff_info_array      = split(";",$gff_info);
        my $gff_ID              = $gff_info_array[0];
        my $gff_XM              = $gff_info_array[2];

        $gff_ID                 =~s/ID=//;
        $gff_XM                 =~s/.+Genbank://;
        if (not exists $gff_hash{$gff_ID}){
                $gff_hash{$gff_ID}      = $gff_XM;
        }
        else{
                print STDERR "GFF ID already in use: $gff_ID - $gff_XM\n";
        }
}
close(GFF) || die;



my %gene_hash;
open(CUFF,"<",$cuff_file) || die;
my $cuff_header = <CUFF>;
while(<CUFF>){
        chomp;
        my $cuff_line   = $_;
        my @cuff_split  = split("\t",$cuff_line);
        my $rna_ID      = $cuff_split[0];
        my $gene_ID     = $cuff_split[3];
        my $length      = $cuff_split[7];
        my $egg_fpkm    = $cuff_split[9];
        my $five_fpkm   = $cuff_split[13];
        my $lli_fpkm    = $cuff_split[17];
        next unless (exists $gff_hash{$rna_ID});
        my $mrna_ID     = $gff_hash{$rna_ID};
	print "$egg_fpkm,1,$mrna_ID\n";
	print "$five_fpkm,4,$mrna_ID\n";
	print "$lli_fpkm,5,$mrna_ID\n";
}
close(CUFF) ||die;



