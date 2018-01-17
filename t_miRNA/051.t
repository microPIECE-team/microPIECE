#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
#nearly all parts of the script are simple pipeline calls 
#still it is to check if the final computation of RPM values is correct.
mkdir -p 050/t_con1_out/
my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/051_overview_scan.pl',['050/t_con1/ 050/t_DB/t_genome.fa 050/t_DB/t_mirna_precursor.fa 050/t_DB/t_ncRNA.fa 050/t_DB/t_mrna.fa 050/t_con1_out/ 100 bwa 4']);


open(HP,"<","050/t_con1_out/mapping.csv") || die;


my @got_mature	= <MA>;
my @got_hairpin	= <HP>;

#my $got = &parser($stdout);
#my $expected = &parser(join('',<DATA>));
my @exp_hairpin;
my @exp_mature;

while(<DATA>){
	if(s/^hp\t//){
		push(@exp_hairpin,$_);
	}
	elsif(s/^ma\t//){
		push(@exp_mature,$_);
	}
}



is_deeply(\@got_mature,\@exp_mature,'051 mature output as expected');
is_deeply(\@got_hairpin,\@exp_hairpin,'051 hairpin output as expected');


done_testing();



sub parser{
        my $p_stdout            =       $_[0];
        my $p_header            =       "";
        my %p_hash;
        foreach my $p_line ( split("\n",$p_stdout)){
                if($p_line=~/^>/){
                        $p_header       = $p_line;
                }
                else{
                        $p_hash{$p_header}=$p_line;
                }
        }

}

__DATA__
length;rpm;type
30;1;miRNA
30;10;ncRNA
30;20;mRNA
30;5;other
