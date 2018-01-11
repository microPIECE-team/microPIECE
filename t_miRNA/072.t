#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
use Switch;
# result file ends with an empty line

#my %exp_str	= (
#">tca-mir-9e"=>">tca-mir-9e(-37.00)[tca-miR-9e-5p:27-47][tca-miR-9e-3p:64-85]aacag-caauUGUGguguugagucggauuaucuCUUUGUGACUAGUUUAUgauguu|||||||||||||||||||||||||||||||||||||||ucagcucagcuuaguagAGAAACACUGGUCGAAUAcuauaaccgaguuccuCGUA-u",
#">tca-mir-304"=>">tca-mir-304(-19.70)[tca-miR-304:13-35]augcaauagaAAA-UcuuucuAUCUCAAUUGUAAAGUGUAcuaguua|||||||||||||||||||||||||||agguagaguuaauauuucgcaugauuaaa-ccuuuaccagcaacauuu",
#">tca-mir-3905"=>">tca-mir-3905(-52.50)[tca-miR-3905-5p:32-50]uuaucucuca-ggaa--AAACC----ggguaacgcgcuggcggGGGGGGGGCGCUUGccuucggcc|||||||||||||||||||||||||||||||||||||||gcccauugcgugaccgccccccuuccgcgaacggaagccggaacuc-cucaagggcaaaaaa-ccac");


my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/072_create_mirbase_struct.pl',[qw(070/072_test_hairpin.fa 070/072_test_mature.fa 070/072_miRNA.str)]);

my %got_str	= %{&parser("custom.str")};
my %exp_str	= %{&parser("070/072_reference_miRNA.str")};

#foreach(keys %got_str){
#	print "################\n$_ :\n$got_str{$_}\n";
#}


#foreach(keys %exp_str){
#	if (not exists $got_str{$_}){
#		print "$_ not in GOT hash\n";
#		die;
#	}
#	print "****************\n$_ :\n$got_str{$_}$exp_str{$_}##################\n";
#}

#foreach(keys %exp_str){
#	is($got_str{$_},$exp_str{$_},'072 '.$_);
#}
is_deeply(\%got_str,\%exp_str,'072 output as expected');
done_testing();


sub parser{
        my $p_dat	=       $_[0];
	my %p_dict;
	my $p_header;
	open(DAT,"<",$p_dat) || die;	
	while(<DAT>){
		chomp;
		my $p_line	= $_;
        #foreach my $p_line ( split("\n",$p_dat)){
		if($p_line=~/^>/){
			$p_header= (split(" ",$p_line))[0];
			$p_dict{$p_header}	= "$p_line\n";
		}
		else{
			$p_dict{$p_header}	.= "$p_line\n";
		}		
        }
	close(DAT) || die;
        return (\%p_dict);
}



