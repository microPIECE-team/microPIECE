#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp qw/ :POSIX /;
use Digest::MD5;

my @testcases = (
    {
	binsize => 1,
	bam     => "t/testset_binning/hs10_sim.bam",
	expected => "be9d43db01c40a03dcc5d1e2322e5506",
    },
    {
	binsize => 5,
	bam     => "t/testset_binning/hs10_sim.bam",
	expected => "68e1d39df9c684c47a40dce47ae081df",
    },
    {
	binsize => 20,
	bam     => "t/testset_binning/hs10_sim.bam",
	expected => "980eb941363d095a2085b5d2d93ef7ae",
    },
    {
	binsize => 100,
	bam     => "t/testset_binning/hs10_sim.bam",
	expected => "a09e9034e1cc1e83cd0d39e928a421e9",
    },
    {
	binsize => 250,
	bam     => "t/testset_binning/hs10_sim.bam",
	expected => "96fda737eb1d8211194e1e4bcb691f49",
    },
    {
	binsize => 500,
	bam     => "t/testset_binning/hs10_sim.bam",
	expected => "83f250501e2f603b5b0736cb5721ea0d",
    },
    {
	binsize => 1,
	bam     => "t/testset_binning/hs10_sim.bam,t/testset_binning/hs10_sim.bam",
	expected => "be9d43db01c40a03dcc5d1e2322e5506",
#	expected => "620ce8de27bc861d88434517be21cda5",
    },
    {
	binsize => 5,
	bam     => "t/testset_binning/hs10_sim.bam,t/testset_binning/hs10_sim.bam",
	expected => "68e1d39df9c684c47a40dce47ae081df",
#	expected => "078b8c5d6275574af94f015eb2bc4ab3",
    },
    {
	binsize => 20,
	bam     => "t/testset_binning/hs10_sim.bam,t/testset_binning/hs10_sim.bam",
	expected => "980eb941363d095a2085b5d2d93ef7ae",
#	expected => "f06e16f9075b2950cfabc196b7c11be9",
    },
    {
	binsize => 100,
	bam     => "t/testset_binning/hs10_sim.bam,t/testset_binning/hs10_sim.bam",
	expected => "a09e9034e1cc1e83cd0d39e928a421e9",
#	expected => "e3610027b051eba616c844f0f62f8586",
    },
    {
	binsize => 250,
	bam     => "t/testset_binning/hs10_sim.bam,t/testset_binning/hs10_sim.bam",
	expected => "96fda737eb1d8211194e1e4bcb691f49",
#	expected => "5cecb7cc4948c1c24bb42ce0976c7ae8",
    },
    {
	binsize => 500,
	bam     => "t/testset_binning/hs10_sim.bam,t/testset_binning/hs10_sim.bam",
	expected => "83f250501e2f603b5b0736cb5721ea0d",
#	expected => "6b7b623a998e2778503c9cd9c513f2f7",
    },
    {
	binsize => 1,
	bam     => "t/testset_binning/hs10_sim.bam",
	transcript => "t/testset_binning/tca_testset_small.gff",
	expected => "be9d43db01c40a03dcc5d1e2322e5506",
    },
    {
	binsize => 5,
	bam     => "t/testset_binning/hs10_sim.bam",
	transcript => "t/testset_binning/tca_testset_small.gff",
	expected => "68e1d39df9c684c47a40dce47ae081df",
    },
    {
	binsize => 20,
	bam     => "t/testset_binning/hs10_sim.bam",
	transcript => "t/testset_binning/tca_testset_small.gff",
	expected => "980eb941363d095a2085b5d2d93ef7ae",
    },
    {
	binsize => 100,
	bam     => "t/testset_binning/hs10_sim.bam",
	transcript => "t/testset_binning/tca_testset_small.gff",
	expected => "a09e9034e1cc1e83cd0d39e928a421e9",
    },
    {
	binsize => 250,
	bam     => "t/testset_binning/hs10_sim.bam",
	transcript => "t/testset_binning/tca_testset_small.gff",
	expected => "96fda737eb1d8211194e1e4bcb691f49",
    },
    {
	binsize => 500,
	bam     => "t/testset_binning/hs10_sim.bam",
	transcript => "t/testset_binning/tca_testset_small.gff",
	expected => "83f250501e2f603b5b0736cb5721ea0d",
    },
    );

my $script = 'scripts/CLIP_binned_bed_from_bam_and_transcripts_for_piranha.pl';

for(my $i=0; $i<@testcases; $i++)
{
    my $current_test = $testcases[$i];

    my @param = (
	"--bam", $current_test->{bam},
	"--size", $current_test->{binsize}
	);

    if (exists $current_test->{transcript})
    {
	push(@param, ("--transcripts", $current_test->{transcript}));
    }

    my ($return,$stdout,$stderr)=run_script($script, \@param );
    is(Test::Script::Run::last_script_exit_code(), 0, 'Test for testset '.$i.' should run and return 0');
    #diag(Dumper({param => \@param, stdout => $stdout, stderr => $stderr})); use Data::Dumper;

    my $got	= parser(\$stdout);
    is($got,$current_test->{expected}, "Output for testset ".$i." as expected");
}

done_testing();

sub parser
{
    my ($out) = @_;

    my $ctx = Digest::MD5->new;
    $ctx->add(${$out});
    return $ctx->hexdigest();
}
