#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;


my $expected_data	= qq{>tca-miR-3905-5p
GGAAGGGGGAGCCGCCUUG
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
};



my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/021_parse_miRDeep2_output.pl',["-mirdeep_out", "020/021_test_mirdeep2.csv", "-mature_fasta", "020/021_test_mature.fa", "-precursor_copies", "tca-mir-3811c-1,tca-mir-3811c-2,tca-mir-3851a-1,tca-mir-3851a-2"]);


isnt($stdout,"",'021 output is not empty');
my $got = join("\n",sort(split("\n",$stdout)));
my $expected = join("\n",sort(split("\n",$expected_data)));
is($got,$expected,'021 output as expected');

done_testing();




