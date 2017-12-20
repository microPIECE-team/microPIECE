#! /usr/bin/perl
use strict;
use warnings;
use RNA::HairpinFigure qw/draw/;

#>cel-let-7 (-42.90)   [cel-let-7-5p:17-38] [cel-let-7-3p:60-81]
#
#------uaca    gga             U              ---  aaua 
#          cugu   uccggUGAGGUAG AGGUUGUAUAGUUu   gg    u
#          ||||   ||||||||||||| ||||||||||||||   ||     
#          gaca   aggCCAUUCCAUC UUUAACGUAUCaag   cc    u
#agcuucucaa    --g             U              ugg  acca 


#my $name   = 'hsa-mir-92a-1 MI0000093 Homo sapiens miR-92a-1 stem-loop';
#my $seq    = 'CUUUCUACACAGGUUGGGAUCGGUUGCAAUGCUGUGUUUCUGUAUGGUAUUGCACUUGUCCCGGCCUGUUGAGUUUGG';
#my $struct = '..(((...((((((((((((.(((.(((((((((((......)))))))))))))).)))))))))))).))).....';
#my $figure = draw( $seq, $struct );
#print ">$name\n$seq\n$struct\n$figure\n";


#>NC_007417.3|344915..344975|-|126.9|result_bwt1.csvresult_bwa.csv|ame-miR-6000a-3p|high_conf
#GGUUGGCAUAAGGUGGUACCAUGUAACAUUUUAACCCAUAGUACGACCCAUGCCGACUCA
#(((((((((..(((.((((.(((.............))).)))).))).))))))))).. (-21.92)



my $hairpin_file	= $ARGV[0];
my $mature_file		= $ARGV[1];
my $rnafold		= "RNAfold -noPS";


my %hairpin_hash	= %{&read_fasta($hairpin_file)};
my %mature_hash		= %{&read_fasta($mature_file)};


open(OUT,">","custom.str") || die;

foreach(keys %hairpin_hash){
	my $hairpin_id		= $_;
	my $hairpin_seq		= $hairpin_hash{$hairpin_id};
	if(not exists $mature_hash{"$hairpin_id-5p"}){
		print STDERR "check $hairpin_id-5p manually ! Did not find hairpin!\n";
	}
	elsif(not exists $mature_hash{"$hairpin_id-3p"}){
		print STDERR "check $hairpin_id-3p manually ! Did not find hairpin\n";
	}
	else{
		my $mature5p_id		= "$hairpin_id-5p";
		$mature5p_id		=~s/>//;
		my $mature3p_id		= "$hairpin_id-3p";
		$mature3p_id		=~s/>//;
		my $mature5p_seq	= $mature_hash{">$mature5p_id"};
		my $mature3p_seq	= $mature_hash{">$mature3p_id"};
		my $mature5p_start	= index($hairpin_seq,$mature5p_seq)+1;
		my $mature5p_stop	= $mature5p_start + length($mature5p_seq)-1;
		my $mature3p_start	= index($hairpin_seq,$mature3p_seq)+1;
		my $mature3p_stop	= $mature3p_start + length($mature3p_seq)-1;
		
		open(TMP,">","tmp_hairpin.fa") || die;
		print TMP "$hairpin_id\n$hairpin_seq";
		close(TMP) || die;
		system("$rnafold < tmp_hairpin.fa > tmp_struct.rna");
		my $hairpin_struct;	#((..))
		my $hairpin_energy;	#(-20.00)
		open(RNA,"<","tmp_struct.rna") || die;
		while(<RNA>){
			chomp;
			my $rna_line 	= $_;
			next if ($.==1);
			next if ($.==2);
			my @rna_split	= split(" ",$rna_line);
			$hairpin_struct	= $rna_split[0];
			$hairpin_energy	= $rna_split[1];
		}
		close(RNA) || die;
		my $hairpin_figure	= draw($hairpin_seq,$hairpin_struct);
		print OUT "$hairpin_id $hairpin_energy   [$mature5p_id:$mature5p_start-$mature5p_stop] [$mature3p_id:$mature3p_start-$mature3p_stop]";
		print OUT "\n\n";
		print OUT "$hairpin_figure\n\n";
	}	
}

close(OUT) || die;
 









# {mirID}       = seq;
sub read_fasta{
        my $rf_file     = $_[0];
        my %rf_hash;    #{mirID} = seq	
	my $rf_header	= "";
        open(RF,"<",$rf_file) || die;
        while(<RF>){
                chomp;
		my $rf_line             = $_;
                my @rf_split            = split(" ",$rf_line);
                if(/^>/){
			$rf_header	= $rf_split[0];
			$rf_header	=lc($rf_header);
			$rf_hash{$rf_header}	= "";
		}
		else{
			$rf_hash{$rf_header}	.= $rf_line;
		}
        }
        close(RF) || die;
        return(\%rf_hash);
}

