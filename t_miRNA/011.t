#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;



my ($return,$stdout,$stderr)=run_script('../scripts/miRNA/011_mirbase_files.pl',["-species", "tca", "-precursor_file", "010/hairpin.fa", "-mature_file", "010/mature.fa", "-organism", "010/organisms.txt", "-out", "010/"]);



####################################

my $got_data    = &parser("010/tca_mature_mirbase.fa");
my $exp_data    = &parser("010/011_tca_mature_mirbase.fa");

my $got = join("\n",sort(split("\n",$got_data)));
my $exp = join("\n",sort(split("\n",$exp_data)));
is($got,$exp,'011 species mature output as expected');

####################################
####################################

$got_data    = &parser("010/tca_precursor_mirbase.fa");
$exp_data    = &parser("010/011_tca_precursor_mirbase.fa");

$got = join("\n",sort(split("\n",$got_data)));
$exp = join("\n",sort(split("\n",$exp_data)));
is($got,$exp,'011 species precursor output as expected');

####################################
####################################

$got_data    = &parser("010/mature.fa-no-tca.fa");
$exp_data    = &parser("010/011_mature.fa-no-tca.fa");

$got = join("\n",sort(split("\n",$got_data)));
$exp = join("\n",sort(split("\n",$exp_data)));
is($got,$exp,'011 other species mature output as expected');

####################################

done_testing();


sub parser{
        my $p_dat       =       $_[0];
        my $p_return    = "";
        my $p_header;
        open(DAT,"<",$p_dat) || die;
        local $/;
        $p_return       = <DAT>;
        close(DAT) || die;
        return ($p_return);
}

