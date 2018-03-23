#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp qw/ :POSIX /;
use Digest::MD5;

my @testcases = (
    {
	input => "MINING_curate_mirdeep2fasta_case1.dat",
	expected => "b5f303bb04b2a8947652bcab1a9171f4",
    },
    {
	input => "MINING_curate_mirdeep2fasta_case2.dat",
	expected => "b5f303bb04b2a8947652bcab1a9171f4",
    }
    );

my $mature  = tmpnam();
my $hairpin = tmpnam();

my ($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl');
isnt(Test::Script::Run::last_script_exit_code(), 0, 'Without arguments is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify --csv file/, "Test for missing csv file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input} ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With file for query is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --cutoff cutoffscore/, "Test for missing cutoff");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10 ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With value for cufoff is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --matureout file/, "Test for missing matureout file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10, "--matureout", $mature ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With value for mature is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --hairpinout file/, "Test for missing hairpinout file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10, "--matureout", $mature, "--hairpinout", $hairpin ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With value for hairpin is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --species three-letter-species-code/, "Test for missing species");


foreach my $current_test (@testcases)
{
    my $mature  = tmpnam();
    my $hairpin = tmpnam();

    my ($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', [
						"--csv", $current_test->{input}, 
						"--cutoff", 10, 
						"--matureout", $mature, 
						"--hairpinout", $hairpin, 
						"--species", "tca"
					    ] );
    is(Test::Script::Run::last_script_exit_code(), 0, 'With value for hairpin is should exit with exit code not equals to 0');
    diag($stdout);
    diag($stderr);

    my $got	= parser($mature, $hairpin);
    is($got,$current_test->{expected}, "output for input: '".$current_test->{input}."' as expected");
}

done_testing();

sub parser
{
    my @files = @_;

    my $ctx = Digest::MD5->new;
    
    foreach my $file (@files)
    {
	open(my $fh, "<", $file) || die "Unable to open input file '$file': $!\n";
	#	$ctx->add($fh);
	while (<$fh>)
	{
	    $ctx->add($_);
	}
	close($fh) || die "Unable to close input file '$file': $!\n";
    }
    
    return $ctx->hexdigest();
}
