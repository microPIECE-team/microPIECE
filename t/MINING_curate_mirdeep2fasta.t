#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp qw/ :POSIX /;
use Digest::MD5;

my @testcases = (
    {
	input => "t/MINING_curate_mirdeep2fasta_case1.dat",
	expected => "a545dbc1e0ea240c00dd54032707ed44",
    },
    {
	input => "t/MINING_curate_mirdeep2fasta_case2.dat",
	expected => "a545dbc1e0ea240c00dd54032707ed44",
    },
    {
	input => "t/MINING_curate_mirdeep2fasta_case3.dat",
	expected => "a545dbc1e0ea240c00dd54032707ed44",
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

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}."non_existing" ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With non existing file for query is should exit with exit code not equals to 0');
like($stderr, qr/Input file does not exist/, "Test for non-existing input file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10 ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With value for cufoff is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --matureout file/, "Test for missing matureout file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10, "--matureout", $mature ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With value for mature is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --hairpinout file/, "Test for missing hairpinout file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10, "--matureout", $testcases[0]{input} ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With an existing file for mature sequences is should exit with exit code not equals to 0');
like($stderr, qr/Mature file exists and will not be overwritten/, "Test for existing mature file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10, "--matureout", $mature, "--hairpinout", $hairpin ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With value for hairpin is should exit with exit code not equals to 0');
like($stderr, qr/Need to specify --species three-letter-species-code/, "Test for missing species");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", $testcases[0]{input}, "--cutoff", 10, "--matureout", $mature, "--hairpinout", $testcases[0]{input} ] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With an existing file for hairpin sequences is should exit with exit code not equals to 0');
like($stderr, qr/Hairpin file exists and will not be overwritten/, "Test for existing hairpin file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', ["--csv", "t/MINING_curate_mirdeep2fasta_collision.dat", "--cutoff", 10, "--matureout", $mature, "--hairpinout", $hairpin, "--species", "tca" ] );
is(Test::Script::Run::last_script_exit_code(), 0, 'Collision should be indicated, but due to same sequence, it should not cause the exit code not equals to 0');
like($stderr, qr/Collision detected, but assuming to have found a genomic copy/, "Test for collision in numbering");
unlink($mature) || die "Unable to unlink mature file '$mature': $!\n";
unlink($hairpin) || die "Unable to unlink hairpin file '$hairpin': $!\n";

foreach my $current_test (@testcases)
{
    my ($return,$stdout,$stderr)=run_script('../scripts/MINING_curate_mirdeep2fasta.pl', [
						"--csv", $current_test->{input}, 
						"--cutoff", 10, 
						"--matureout", $mature, 
						"--hairpinout", $hairpin, 
						"--species", "tca"
					    ] );
    is(Test::Script::Run::last_script_exit_code(), 0, 'With input values it should run and return 0');

    my $got	= parser($mature, $hairpin);
    is($got,$current_test->{expected}, "output for input: '".$current_test->{input}."' as expected");
    unlink($mature) || die "Unable to unlink mature file '$mature': $!\n";
    unlink($hairpin) || die "Unable to unlink hairpin file '$hairpin': $!\n";
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
