#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
# result file ends with an empty line


my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/072_create_mirbase_struct.pl',["-hairpin","070/072_test_hairpin.fa","-mature", "070/072_test_mature.fa","-struct", "070/072_miRNA.str"]);

my $got_data	= &parser("custom.str");
my $exp_data	= &parser("070/072_reference_miRNA.str");

my $got = join("\n",sort(split("\n",$got_data)));
my $exp = join("\n",sort(split("\n",$exp_data)));



is($got,$exp,'072 output as expected');
done_testing();


sub parser{
        my $p_dat	=       $_[0];
	my $p_return	= "";
	my $p_header;
	open(DAT,"<",$p_dat) || die;	
	while(<DAT>){
		chomp;
		my $p_line	= $_;
		$p_return	.= "$p_line\n";
        }
	close(DAT) || die;
        return ($p_return);
}



