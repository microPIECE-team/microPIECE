#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

my $first_file = qq{seq	name	freq	mir	start	end	mism	add	t5	t3	s5	s3	DB	precursor	ambiguity
};

my $second_file = qq{seq	name	freq	mir	start	end	mism	add	t5	t3	s5	s3	DB	precursor	ambiguity
};

my $expected_output = qq{
};

my ($fh1, $filename1) = File::Temp::tempfile();

print $fh1 $first_file;
close($fh1) || die;

my ($fh2, $filename2) = File::Temp::tempfile();

print $fh2 $second_file;
close($fh2) || die;

my (undef, $filename3) = File::Temp::tempfile(OPEN => 0);

my ($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl');
isnt(Test::Script::Run::last_script_exit_code(), 0, 'Without arguments is should exit with exit code not equals to 0';

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl', ["--query", $filename1] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With file for query is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a subject FASTA file via --subject parameter/, "Test for missing subject file");

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl', ["--query", $filename1."not_there"] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With non existing file for query is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a query FASTA file via --query parameter/, "Test with non-existing query file");

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl', ["--subject", $filename2] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With file for subject is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a query FASTA file via --query parameter/, "Test for missing query file");

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl', ["--query", $filename1, "--subject", $filename2."not_there"] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With non existing subject for query is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a subject FASTA file via --subject parameter/, "Test with non-existing missing subject file");

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl',["--query", $filename1, "--subject", $filename2]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With files for query and subject it should exit with exit code 0');

my $got	= parser($stdout);
my $expected = parser($expected_output);
is_deeply($got,$expected,'output as expected');

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl',["--query", $filename1, "--subject", $filename2, "--out", $filename3]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With files for query, subject, and out it should exit with exit code 0');

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl',["--query", $filename1, "--subject", $filename2, "--out", $filename3]);
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With files for query, subject, and existing out it should exit without exit code 0');

($return,$stdout,$stderr)=run_script('../scripts/ISOMIR_reformat_isomirs.pl',["--query", $filename1, "--subject", $filename2, "--out", $filename3, "--overwrite"]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With files for query, subject, and existing out and overwrite it should exit with exit code 0');

done_testing();

sub parser
{
    my ($in) = @_;

    # split by newline and ignore empty lines
    my @lines = grep {$_ !~ /^\s*$/ } (split(/\n/, $in));

    # split each line on white spaces
    for(my $i=0; $i<@lines; $i++)
    {
	$lines[$i] = [split(/\s+/, $lines[$i])];
    }

    # sort the array at fields 0&1 (cmp), 6,7,8,9 (<=>)
    @lines = sort { $a->[0] cmp $b->[0]
			||
			$a->[1] cmp $b->[1]
			||
			$a->[6] <=> $b->[6]
			||
			$a->[7] <=> $b->[7]
			||
			$a->[8] <=> $b->[8]
			||
			$a->[9] <=> $b->[9]
    } (@lines);

    return \@lines;
}
