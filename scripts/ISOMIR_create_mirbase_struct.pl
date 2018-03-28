#!/usr/bin/env perl
use strict;
use warnings;
use RNA::HairpinFigure qw/draw/;
use Getopt::Long;
use File::Temp qw(tmpnam);

use FindBin;
use lib "$FindBin::Bin/lib";

use mining;

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

my $output;
my $mirbase_dat;

GetOptions(
    "mirbasedat=s" => \$mirbase_dat,
    "out=s"        => \$output
    ) || die;

die "Need to specify a existing mirbase file via --mirbasedat parameter\n" unless (defined $mirbase_dat && -e $mirbase_dat);
die "Need to specify a non-existing output file via --out parameter\n" unless (defined $output && ( ! -e $output ));

my $rnafold		= "RNAfold -noPS";

my $mirnas = mining::parse_mirbase_dat($mirbase_dat, undef);

my $foldinput = tmpnam();

open(OUT,">",$output) || die "Unable to open file '$output' for writing: $!\n";

foreach my $entry (@{$mirnas})
{
    open(FH, ">", $foldinput) || die "Unable to open file '$foldinput': $!\n";
    print FH ">temp\n", $entry->{seq}, "\n";
    close(FH) || die "Unable to close file '$foldinput': $!\n";

    my $cmd = "$rnafold <$foldinput";
    open(RNAFOLD, "$cmd|") || die "Unable to open pipe via '$cmd': $!\n";
    <RNAFOLD>;
    <RNAFOLD>;
    my ($struct, $energy) = split(/\s+/, scalar <RNAFOLD>);
    close(RNAFOLD) || die "Unable to close pipe: $!\n";

    # generate the plot
    my $seq = lc($entry->{seq});
    my @mature_infos = ();
    foreach my $mature (@{$entry->{matures}})
    {
	my $start = $mature->{start}-1;
	my $len   = $mature->{stop}-$mature->{start}+1;
	substr($seq, $start, $len, uc(substr($seq, $start, $len)));

	push(@mature_infos, sprintf("[%s:%d-%d]", $mature->{name}, $mature->{start}, $mature->{stop}));

    }

    printf OUT ">%s %s %s\n\n", $entry->{precursor}, $energy, join(" ", @mature_infos);
    print  OUT draw($seq, $struct), "\n\n";
}

close(OUT) || die "Unable to close file '$output' after writing: $!\n";
