#! /usr/bin/perl
use strict;
use warnings;
use Switch;
# loop through each replicate
# count reads by length that mapped to the different targets (genome, ncRNA, miRNA,mRNA...)
# calculate RPM for each replicate target length
# sum up each replicate target length RPM
# divide by 4 (replicates in this case)


# bowtie2-build reference.fastq IDX
# needed for genome, precursor, ncrna and mrna files

my $read_folder         = $ARGV[0];     # per condition with replicates
my $genome_idx          = $ARGV[1];
my $mir_idx 		= $ARGV[2];
my $ncrna_idx           = $ARGV[3];
my $mrna_idx            = $ARGV[4];
my $out_folder          = $ARGV[5];
my $threads             = $ARGV[6];
my $aligner		= $ARGV[7];	#bwt1 bwt2 bwa
my $replicates		= $ARGV[8]; 	# 4 

switch($aligner){
	case	"bwt2"	{print "using BOWTIE2\n";}
	case	"bwt1"	{print "using BOWTIE1\n";}
	case 	"bwa"	{print "using BWA\n";}
	else		{print "No valid aligner! { bwt1 | bwt2 | bwa } Terminating ...\n"; exit;}
}


                                #            unmappable      , miRNA        , ncRNA      , mRNA        , other
my %result_hash;                # {len} = \@(rpm__unal_genome,rpm_aln_mir, rpm_aln_ncRNA, rpm_aln_mrna, rpm_unal_mrna)
# init result_hash
my @range_array         = (17..40);
foreach(@range_array){
        my $range_len   = $_;
        my @dummy_array = (0,0,0,0,0);
        $result_hash{$range_len}=\@dummy_array;
}


# unmappable
# mirna
# ncRNA
# mRNA
# other

my $bowtie2		= "bowtie2 --very-sensitive-local -p $threads"; # -x BWT2_index
my $bowtie1		= "bowtie --best --strata -p $threads";
my $bwa_aln		= "bwa aln -n 1 -o 0 -e 0 -k 1 -t $threads";
my $bwa_samse		= "bwa samse";
my $samtools		= "samtools";
my $bedtools		= "bedtools";
opendir DIR, $read_folder;
my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;       # each fastq file in array without path
closedir DIR;




foreach(@files){
        my $file        = $_;
        my $in_path     = $read_folder.$file;
        my $out_path    = $out_folder.$file;
        my $read_count  = &countREADS($in_path);
        print "mapping to genome\n";
        my $genome_ual  = $out_path."-genome_ual.fq";
        my $genome_aln  = $out_path."-genome_aln.fq";
	my $genome_sam	= $out_path."-genome.sam";
	switch($aligner){

		case	"bwt2"	{
			print "$bowtie2 --un $genome_ual --al $genome_aln -x $genome_idx $in_path -S $genome_sam\n";
			system("$bowtie2 --un $genome_ual --al $genome_aln -x $genome_idx $in_path -S $genome_sam");
		}
		case	"bwt1"	{
			print "$bowtie1 --un $genome_ual -- al $genome_aln $genome_idx $in_path -S $genome_sam\n";
			system("$bowtie1 --un $genome_ual -- al $genome_aln $genome_idx $in_path -S $genome_sam");
		}
		case 	"bwa"	{
			my $genome_sai		= $out_path."-genome.sai";
			my $genome_ual_bam	= $out_path."-genome_ual.bam";
			my $genome_aln_bam	= $out_path."-genome_aln.bam";
			print "$bwa_aln -f $genome_sai $genome_idx $in_path\n";
			system("$bwa_aln -f $genome_sai $genome_idx $in_path");
			print "$bwa_samse -f $genome_sam $genome_idx $genome_sai $in_path\n";
			system("$bwa_samse -f $genome_sam $genome_idx $genome_sai $in_path");
			system("$samtools view -b -f 4 $genome_sam > $genome_ual_bam");
			system("$bedtools bamtofastq -i $genome_ual_bam -fq $genome_ual");
			system("$samtools view -b -F 4 $genome_sam > $genome_aln_bam");
			system("$bedtools bamtofastq -i $genome_aln_bam -fq $genome_aln");
		}
	
	}
        my %genome_hash = %{&countFASTA($genome_ual,$read_count)};
        foreach(keys %genome_hash){
                my $genome_len                  = $_;
                my $genome_rpm                  = $genome_hash{$genome_len};
                my @genome_array                = @{$result_hash{$genome_len}};
                $genome_array[0]               += $genome_rpm;
                $result_hash{$genome_len}       = \@genome_array;
        }
        # use genome_aln for next test
        print "mapping ncRNA\n";
        my $ncrna_ual   = $out_path."-ncrna_ual.fq";
        my $ncrna_aln   = $out_path."-ncrna_aln.fq";
	my $ncrna_sam	= $out_path."-ncrna.sam";
	switch($aligner){
		case	"bwt2"	{
			print "$bowtie2  --un $ncrna_ual --al $ncrna_aln -x $ncrna_idx $genome_aln -S $ncrna_sam\n";
			system("$bowtie2  --un $ncrna_ual --al $ncrna_aln -x $ncrna_idx $genome_aln -S $ncrna_sam");
		}
		case	"bwt1"	{
                        print "$bowtie1 --un $ncrna_ual -- al $ncrna_aln $ncrna_idx $genome_aln -S $ncrna_sam\n";
                        system("$bowtie1 --un $ncrna_ual -- al $ncrna_aln $ncrna_idx $genome_aln -S $ncrna_sam");
		}
		case	"bwa"	{
                        my $ncrna_sai  		= $out_path."-ncrna.sai";
			my $ncrna_ual_bam	= $out_path."-ncrna_ual.bam";
			my $ncrna_aln_bam	= $out_path."-ncrna_aln.bam";
                        print "$bwa_aln -f $ncrna_sai $ncrna_idx $genome_aln\n";
                        system("$bwa_aln -f $ncrna_sai $ncrna_idx $genome_aln");
                        print "$bwa_samse -f $ncrna_sam $ncrna_idx $ncrna_sai $genome_aln\n";
                        system("$bwa_samse -f $ncrna_sam $ncrna_idx $ncrna_sai $genome_aln");
			system("$samtools view -b -f 4 $ncrna_sam > $ncrna_ual_bam");
			system("$bedtools bamtofastq -i $ncrna_ual_bam -fq $ncrna_ual");
			system("$samtools view -b -F 4 $ncrna_sam > $ncrna_aln_bam");
			system("$bedtools bamtofastq -i $ncrna_aln_bam -fq $ncrna_aln");
		}
	}
        my %ncrna_hash  = %{&countFASTA($ncrna_aln,$read_count)};
        foreach(keys %ncrna_hash){
                my $ncrna_len                   = $_;
                my $ncrna_rpm                   = $ncrna_hash{$ncrna_len};
                my @ncrna_array                 = @{$result_hash{$ncrna_len}};
                $ncrna_array[2]                += $ncrna_rpm;
                $result_hash{$ncrna_len}        = \@ncrna_array;
        }


        # --norc
        # use ncrna_ual for next test
        print "mapping miRNA\n";
        my $mir_ual     = $out_path."-mir_ual.fq";
        my $mir_aln     = $out_path."-mir_aln.fq";
	my $mir_sam	= $out_path."-mir.sam";

	switch($aligner){
		case	"bwt2"	{
			print "$bowtie2  --un $mir_ual --al $mir_aln -x $mir_idx $ncrna_ual -S $mir_sam\n";
			system("$bowtie2  --un $mir_ual --al $mir_aln -x $mir_idx $ncrna_ual -S $mir_sam");
		}
		case	"bwt1"	{
                        print "$bowtie1 --un $mir_ual -- al $mir_aln $mir_idx $ncrna_ual -S $mir_sam\n";
                        system("$bowtie1 --un $mir_ual -- al $mir_aln $mir_idx $ncrna_ual -S $mir_sam");			
		}
		case	"bwa"	{
                        my $mir_sai  	= $out_path."-mir.sai";
			my $mir_ual_bam	= $out_path."-mir_ual.bam";
			my $mir_aln_bam	= $out_path."-mir_aln.bam";
                        print "$bwa_aln -f $mir_sai $mir_idx $ncrna_ual\n";
                        system("$bwa_aln -f $mir_sai $mir_idx $ncrna_ual");
                        print "$bwa_samse -f $mir_sam $mir_idx $mir_sai $ncrna_ual\n";
                        system("$bwa_samse -f $mir_sam $mir_idx $mir_sai $ncrna_ual");
			system("$samtools view -b -f 4 $mir_sam > $mir_ual_bam");
			system("$bedtools bamtofastq -i $mir_ual_bam -fq $mir_ual");
			system("$samtools view -b -F 4 $mir_sam > $mir_aln_bam");
			system("$bedtools bamtofastq -i $mir_aln_bam -fq $mir_aln");
		}
	}
        my %mir_hash    = %{&countFASTA($mir_aln,$read_count)};
        foreach(keys %mir_hash){
                my $mir_len                     = $_;
                my $mir_rpm                     = $mir_hash{$mir_len};
                my @mir_array                   = @{$result_hash{$mir_len}};
                $mir_array[1]                  += $mir_rpm;
                $result_hash{$mir_len}          = \@mir_array;
        }
        # use mirna_ual for next test
        print "mapping mRNA and others\n";
        my $mrna_ual    = $out_path."-mrna_ual.fq";
        my $mrna_aln    = $out_path."-mrna_aln.fq";
	my $mrna_sam	= $out_path."-mrna.sam";

	switch($aligner){
		case "bwt2"	{
			print "$bowtie2 --un $mrna_ual --al $mrna_aln -x $mrna_idx $mir_ual -S $mrna_sam\n";
			system("$bowtie2 --un $mrna_ual --al $mrna_aln -x $mrna_idx $mir_ual -S $mrna_sam");
		}
		case "bwt1"	{
                        print "$bowtie1 --un $mrna_ual -- al $mrna_aln $mrna_idx $mir_ual -S $mrna_sam\n";
                        system("$bowtie1 --un $mrna_ual -- al $mrna_aln $mrna_idx $mir_ual -S $mrna_sam");
		}
		case "bwa"	{
                        my $mrna_sai  		= $out_path."-mrna.sai";
			my $mrna_ual_bam	= $out_path."-mrna_ual.bam";
			my $mrna_aln_bam	= $out_path."-mrna-aln.bam";
                        print "$bwa_aln -f $mrna_sai $mrna_idx $mir_ual\n";
                        system("$bwa_aln -f $mrna_sai $mrna_idx $mir_ual");
                        print "$bwa_samse -f $mrna_sam $mrna_idx $mrna_sai $mir_ual\n";
                        system("$bwa_samse -f $mrna_sam $mrna_idx $mrna_sai $mir_ual");
			system("$samtools view -b -f 4 $mrna_sam > $mrna_ual_bam");
			system("$bedtools bamtofastq -i $mrna_ual_bam -fq $mrna_ual");
			system("$samtools view -b -F 4 $mrna_sam > $mrna_aln_bam");
			system("$bedtools bamtofastq -i $mrna_aln_bam -fq $mrna_aln");
		}	
	}
        my %mrna_hash   = %{&countFASTA($mrna_aln,$read_count)};
        foreach(keys %mrna_hash){
                my $mrna_len                    = $_;
                my $mrna_rpm                    = $mrna_hash{$mrna_len};
                my @mrna_array                  = @{$result_hash{$mrna_len}};
                $mrna_array[3]                 += $mrna_rpm;
                $result_hash{$mrna_len}         = \@mrna_array;
        }
        # use mrna
        my %other_hash  = %{&countFASTA($mrna_ual,$read_count)};
        foreach(keys %other_hash){
                my $other_len                   = $_;
                my $other_rpm                   = $other_hash{$other_len};
                my @other_array                 = @{$result_hash{$other_len}};
                $other_array[4]                += $other_rpm;
                $result_hash{$other_len}        = \@other_array;
        }
}

my $outfile     = $out_folder."mapping.csv";
open(OUT,">",$outfile);
print OUT "length;rpm;type\n";
foreach(keys %result_hash){
        my $rh_len      = $_;
        my @rh_array    = @{$result_hash{$rh_len}};
        my $rh_genome   = $rh_array[0]/$replicates;
        my $rh_mir      = $rh_array[1]/$replicates;
        my $rh_ncrna    = $rh_array[2]/$replicates;
        my $rh_mrna     = $rh_array[3]/$replicates;
        my $rh_other    = $rh_array[4]/$replicates;
#       print OUT "$rh_len;$rh_genome;'genome'\n";
        print OUT "$rh_len;$rh_mir;miRNA\n";
        print OUT "$rh_len;$rh_ncrna;ncRNA\n";
        print OUT "$rh_len;$rh_mrna;mRNA\n";
        print OUT "$rh_len;$rh_other;other\n";
}
close(OUT);

print "plotting\n";
my $outfilesvg  = $outfile."_plot.svg";
system("./052_stacked_barplot.R $outfile $outfilesvg");







sub countREADS{
        my $cr_file     = $_[0];
        my $cr_file_fa  = $cr_file.".fa";
        system("fastq2fasta.pl $cr_file > $cr_file_fa");
        my $cr_count    = 0;
        open(CR,"<",$cr_file_fa);
        while(<CR>){
                chomp;
                if(/^>/){
                        $cr_count +=1;
                }
        }
        close(CR);
	system("rm -rf $cr_file_fa");
        return($cr_count);
}



# in    : fasta file (alignment output)
# out   : {len} = read count            | {total}       = #all_reads
sub countFASTA{
        my $cf_file     = $_[0];        #fastq
        my $cf_total    = $_[1];        # all reads in lib
        my $cf_file_fa  =$cf_file.".fa";
        system("fastq2fasta.pl $cf_file > $cf_file_fa");
        my %cf_hash;
        open(CF,"<",$cf_file_fa);
        while(<CF>){
                chomp;
                if(not /^>/){
                        my $cf_seq      = $_;
                        my $cf_len      = length($cf_seq);
                        $cf_total       +=1;
                        if(not exists $cf_hash{$cf_len}){
                                $cf_hash{$cf_len}       = 1;
                        }
                        else{
                                $cf_hash{$cf_len}       +=1;
                        }
                }
        }
        close(CF);
        foreach(keys %cf_hash){
                my $cf_key              = $_;
                my $cf_count            = $cf_hash{$cf_key};
                my $cf_rpm              = ($cf_count/$cf_total)*1000000;
                $cf_hash{$cf_key}       = $cf_rpm;
                #print "#### $cf_file ::: len  $cf_key : count $cf_count : total $cf_total | rpm $cf_rpm\n";
        }
        return(\%cf_hash);
}

