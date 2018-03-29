#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp qw/ :POSIX /;
use Digest::MD5;

my @testcases = (
    {
	input => 
q(ID   cel-mir-2         standard; RNA; CEL; 98 BP.
XX
AC   MI0000004;
XX
DE   Caenorhabditis elegans miR-2 stem-loop
XX
RN   [1]
RX   PUBMED; 11679671.
RA   Lau NC, Lim LP, Weinstein EG, Bartel DP;
RT   "An abundant class of tiny RNAs with probable regulatory roles in
RT   Caenorhabditis elegans";
RL   Science. 294:858-862(2001).
XX
RN   [2]
RX   PUBMED; 11679672.
RA   Lee RC, Ambros V;
RT   "An extensive class of small RNAs in Caenorhabditis elegans";
RL   Science. 294:862-864(2001).
XX
RN   [3]
RX   PUBMED; 12672692.
RA   Lim LP, Lau NC, Weinstein EG, Abdelhakim A, Yekta S, Rhoades MW, Burge CB,
RA   Bartel DP;
RT   "The microRNAs of Caenorhabditis elegans";
RL   Genes Dev. 17:991-1008(2003).
XX
RN   [4]
RX   PUBMED; 12747828.
RA   Ambros V, Lee RC, Lavanway A, Williams PT, Jewell D;
RT   "MicroRNAs and other tiny endogenous RNAs in C. elegans";
RL   Curr Biol. 13:807-818(2003).
XX
RN   [5]
RX   PUBMED; 12769849.
RA   Grad Y, Aach J, Hayes GD, Reinhart BJ, Church GM, Ruvkun G, Kim J;
RT   "Computational and experimental identification of C. elegans microRNAs";
RL   Mol Cell. 11:1253-1263(2003).
XX
RN   [6]
RX   PUBMED; 17174894.
RA   Ruby JG, Jan C, Player C, Axtell MJ, Lee W, Nusbaum C, Ge H, Bartel DP;
RT   "Large-scale sequencing reveals 21U-RNAs and additional microRNAs and
RT   endogenous siRNAs in C. elegans";
RL   Cell. 127:1193-1207(2006).
XX
RN   [7]
RX   PUBMED; 19460142.
RA   Kato M, de Lencastre A, Pincus Z, Slack FJ;
RT   "Dynamic expression of small non-coding RNAs, including novel microRNAs
RT   and piRNAs/21U-RNAs, during Caenorhabditis elegans development";
RL   Genome Biol. 10:R54(2009).
XX
RN   [8]
RX   PUBMED; 20062054.
RA   Zisoulis DG, Lovci MT, Wilbert ML, Hutt KR, Liang TY, Pasquinelli AE, Yeo
RA   GW;
RT   "Comprehensive discovery of endogenous Argonaute binding sites in
RT   Caenorhabditis elegans";
RL   Nat Struct Mol Biol. 17:173-179(2010).
XX
RN   [9]
RX   PUBMED; 21307183.
RA   Warf MB, Johnson WE, Bass BL;
RT   "Improved annotation of C. elegans microRNAs by deep sequencing reveals
RT   structures associated with processing by Drosha and Dicer";
RL   RNA. 17:563-577(2011).
XX
DR   RFAM; RF00047; mir-2.
DR   WORMBASE; M04C9/29652-29555; .
XX
FH   Key             Location/Qualifiers
FH
FT   miRNA           20..41
FT                   /accession="MIMAT0020302"
FT                   /product="cel-miR-2-5p"
FT                   /evidence=experimental
FT                   /experiment="Illumina [9]"
FT   miRNA           61..83
FT                   /accession="MIMAT0000004"
FT                   /product="cel-miR-2-3p"
FT                   /evidence=experimental
FT                   /experiment="cloned [1-4], PCR [5], 454 [6], Illumina
FT                   [7,9], CLIPseq [8]"
XX
SQ   Sequence 98 BP; 27 A; 19 C; 22 G; 0 T; 30 other;
     uaaacaguau acagaaagcc aucaaagcgg ugguugaugu guugcaaauu augacuuuca        60
     uaucacagcc agcuuugaug ugcugccugu ugcacugu                                98
//
),
	expected => 
q(>cel-mir-2 (-38.40)   [cel-miR-2-5p:20-41] [cel-miR-2-3p:61-83]

uaa     -au    aa   -         -   GGU       -uu -  a 
   acagu   acag  agc CAUCAAAGC GGU   UGAUGUG   g ca a
   |||||   ||||  ||| ||||||||| |||   |||||||   | || u
   uguca   uguc  uCG GUAGUUUCG CCG   ACUAUac   c gu u
---     cgu    cg   U         A   -AC       uuu a  a 

),
    },
    {
	input =>
q(ID   tca-mir-1    standard; RNA; TCA; 72 BP.
XX
FH   Key             Location/Qualifiers
FH
FT   miRNA           11..32
FT                   /product="tca-miR-1-5p"
FT   miRNA           47..68
FT                   /product="tca-miR-1-3p"
XX
SQ   Sequence 72 BP; 19 A; 14 C; 17 G; 0 T; 22 other;
     uucgcgaguu ccgugcuucc uuacuuccca uagugcuuua aacguaugga auguaaagaa        60
     guauggagcg aa                                                            72
//
),
	expected =>
q(>tca-mir-1 (-28.20)   [tca-miR-1-5p:11-32] [tca-miR-1-3p:47-68]

uucgcga            C    UUC     -  gc
       guuCCGUGCUUC UUAC   CCAUA gu  u
       |||||||||||| ||||   ||||| ||  u
       cGAGGUAUGAAG AAUG   GGUau ca  u
----aag            A    UAA     g  aa

),
    },
    );

my $dat  = tmpnam();
my $out = tmpnam();

open(FH, ">", $dat) || die;
print FH $testcases[0]{input};
close(FH) || die;

my $script = "../scripts/ISOMIR_create_mirbase_struct.pl";

my ($return,$stdout,$stderr)=run_script($script);
isnt(Test::Script::Run::last_script_exit_code(), 0, 'Without arguments is should exit with exit code not equals to 1');
#like($stderr, qr/Need to specify --csv file/, "Test for missing csv file");
foreach ($dat) { unlink($_) || die "Unable to delete '$_': $!\n"; }

foreach my $testcase (@testcases)
{
    open(FH, ">", $dat) || die;
    print FH $testcase->{input};
    close(FH) || die;

    ($return,$stdout,$stderr)=run_script($script, [
					     "--mirbasedat", $dat,
					     "--out", $out
					 ] );
    is(Test::Script::Run::last_script_exit_code(), 0, 'With input values it should run and return 0');
    my $got = "";
    open(FH, "<", $out) || die;
    while (<FH>) { $got .= $_; }
    close(FH) || die;

    my $expected = $testcase->{expected};

    $got      =~ s/\s+/ /g;
    $expected =~ s/\s+/ /g;

    is($got,$expected, "output as expected");
    foreach ($dat, $out) { unlink($_) || die "Unable to delete '$_': $!\n"; }
}

done_testing();
