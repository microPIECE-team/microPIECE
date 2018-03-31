use strict;
use warnings;

use Test::More;

use FindBin;
use lib "$FindBin::Bin/../scripts/lib";

use mining;
use RNA::HairpinFigure;

my @testcases = (
    {
	sequence  => "UUCGCGAGUUCCGUGCUUCCUUACUUCCCAUAGUGCUUUAAACGUAUGGAAUGUAAAGAAGUAUGGAGCGAA",
	structure => ".......((((((((((((.((((...(((((((.......)).)))))...)))).))))))))))))...",
	expected  => "UUCGCGA            C    UUC     -  GC \n".
                     "       GUUCCGUGCUUC UUAC   CCAUA GU  U\n".
                     "       |||||||||||| ||||   ||||| ||  U\n".
                     "       CGAGGUAUGAAG AAUG   GGUAU CA  U\n".
                     "----AAG            A    UAA     G  AA ",
	title     => "Wrong structure from RNA::HairpinFigure"
    },
    {
	sequence  => "uaaacaguauacagaaagcCAUCAAAGCGGUGGUUGAUGUGuugcaaauuaugacuuucaUAUCACAGCCAGCUUUGAUGUGCugccuguugcacugu",
	structure => "...(((((..((((..(((((((((((((((...(((((((..(((.....)).)...)))))))..))).))))))))).)))..))))...)))))",
	expected  => "uaa     -au    aa   -         -   GGU       -uu -  a \n".
                     "   acagu   acag  agc CAUCAAAGC GGU   UGAUGUG   g ca a\n".
                     "   |||||   ||||  ||| ||||||||| |||   |||||||   | || u\n".
                     "   uguca   uguc  uCG GUAGUUUCG CCG   ACUAUac   c gu u\n".
                     "---     cgu    cg   U         A   -AC       uuu a  a ",
	title     => "Correct structure from RNA::HairpinFigure"
    },
);

foreach my $case (@testcases)
{
    my $got = mining::fix_hairpin(draw($case->{sequence}, $case->{structure}), $case->{sequence});
    is($got, $case->{expected}, $case->{title});
}

done_testing();
