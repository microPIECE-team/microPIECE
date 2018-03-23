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
    );

my $script = 'scripts/CLIP_binned_bed_from_bam_and_transcripts_for_piranha.pl';

for(my $i=0; $i<@testcases; $i++)
{
    my $current_test = $testcases[$i];
    
    my ($return,$stdout,$stderr)=run_script($script, [
						"--bam", $current_test->{bam}, 
						"--size", $current_test->{binsize} 
					    ] );
    is(Test::Script::Run::last_script_exit_code(), 0, 'Test for testset '.$i.' should run and return 0');

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
