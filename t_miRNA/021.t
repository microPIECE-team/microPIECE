#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;
use Switch;

my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/021_parse_miRBaseOUT.pl',[qw(020/021_test_mirdeep2.csv 020/021_test_mature.fa tca-mir-3811c-1,tca-mir-3811c-2,tca-mir-3851a-1,tca-mir-3851a-2)]);


my $got = &parser($stdout);
my $expected = &parser(join('',<DATA>));
is_deeply($got,$expected,'021 output as expected');

done_testing();





sub parser{
        my $p_stdout            =       $_[0];
        my $p_header            =       "";
        my %p_hash;
        foreach my $p_line ( split("\n",$p_stdout)){
		if($p_line=~/^>/){
			$p_header	= $p_line;
		}
		else{
			$p_hash{$p_header}=$p_line;
		}
	}

}

__DATA__
>tca-miR-3811c-1-5p
UUCUUGCUCAGGGUUUACAUGU
>tca-miR-3811c-2-5p
UUCUUGCUCAGGGUUUACAUGU
>tca-miR-3811c-1-3p
UAUGUACAACCGUGAGUGGAGACC
>tca-miR-3811c-2-3p
UAUGUACAACCGUGAGUGGAGACC
>tca-miR-3851a-1-3p
UAGUACAUUCCGAAUCGACUCG
>tca-miR-3851a-2-3p
UAGUACAUUCCGAAUCGACUCG
>tca-miR-3851a-1-5p
AGUUGAUUUGGAGUAGUACUAUG
>tca-miR-3851a-2-5p
AGUUGAUUUGGAGUAGUACUAUG
>tca-let-7-3p
CUGUACAGCCUGCUAACUUUCCC
>tca-let-7-5p
UGAGGUAGUAGGUUGUAUAG
>tca-miR-iab-4-5p
ACGUAUACUGAAUGUAUCCUGA
>tca-miR-iab-4-3p
CGGUAUACCUUCAGUAUACGUAAC
>tca-miR-981-5p
CGGGUUUCGGGGCAUCGGAACC
>tca-miR-981-3p
UUCGUUGUCGACGAAACCUGCA
>tca-miR-315-5p
UUUUGAUUGUUGCUCAGAAAGCC
>tca-miR-315-3p
CUUUCGGGCAAUAAUCAUUUCC
