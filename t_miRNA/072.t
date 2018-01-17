#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
# result file ends with an empty line


my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/072_create_mirbase_struct.pl',["-hairpin","070/072_test_hairpin.fa","-mature", "070/072_test_mature.fa","-struct", "070/072_miRNA.str"]);

my %got_str	= %{&parser("custom.str")};
my %exp_str	= %{&parser("070/072_reference_miRNA.str")};



is(%got_str,%exp_str,'072 output as expected');
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



