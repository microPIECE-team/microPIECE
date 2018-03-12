#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

use File::Temp;

my $first_file = qq{>seqA01
CGTTACTT
>seqA02
CGACATACTC
>seqA03
ACTCCCCTGGTGATCTAAAT
>seqA04
GCGACATACTC
>seqA05
ACTCCCCTGGG
};

my $second_file = qq{>seqB01 # 9 perfect (fail)
CGTTACTT
>seqB02 # 10 perfect (pass)
CGACATACTC
>seqB03 # 20 perfect (pass)
ACTCCCCTGGTGATCTAAAT
>seqB04 # 10 perfect, 1 mm, 9 matches (pass)
ACTCCCCTGGAGATCTAAAT
>seqB05 # 10 perfect, 2 mm, 8 matches (fail)
ACTCCCCTGGAAATCTAAAT
>seqB06 # 10 perfect, 1 gap, 9 matches (pass)
ACTCCCCTGGGATCTAAAT
>seqB07 # 10 perfect, 1 gap, 1 match, 1 gap, 7 matches (fail)
ACTCCCCTGGGTCTAAAT
>seqB07 # 20 mm (fail)
GAAGGGCCAAGAAGGGCCAA
>seqB08 # 1 mm, 10 matches (fail)
CCGACATACTC
>seqB09 # 5th mismatch, 19 matches
ACTCTCCTGGTGATCTAAAT
};

my $expected_output = qq{
seqA02	seqB02	100.00	10	0	0	1	10	1	10	0.001	19.6	CGACATACTC	CGACATACTC	10	10	100.0	100.0
seqA03	seqB03	100.00	20	0	0	1	20	1	20	6e-09	38.1	ACTCCCCTGGTGATCTAAAT	ACTCCCCTGGTGATCTAAAT	20	20	100.0	100.0
seqA03	seqB04	95.00	20	1	0	1	20	1	20	3e-07	32.5	ACTCCCCTGGTGATCTAAAT	ACTCCCCTGGAGATCTAAAT	20	20	100.0	100.0
seqA03	seqB06	95.00	20	0	1	1	20	1	19	1e-06	30.7	ACTCCCCTGGTGATCTAAAT	ACTCCCCTGG-GATCTAAAT	20	19	100.0	100.0
seqA05	seqB07	100.00	11	0	0	1	11	1	11	3e-04	21.4	ACTCCCCTGGG	ACTCCCCTGGG	11	18	100.0	61.1
seqA05	seqB06	100.00	11	0	0	1	11	1	11	3e-04	21.4	ACTCCCCTGGG	ACTCCCCTGGG	11	19	100.0	57.9
seqA05	seqB05	100.00	10	0	0	1	10	1	10	0.001	19.6	ACTCCCCTGG	ACTCCCCTGG	11	20	90.9	50.0
seqA05	seqB04	100.00	10	0	0	1	10	1	10	0.001	19.6	ACTCCCCTGG	ACTCCCCTGG	11	20	90.9	50.0
seqA05	seqB03	100.00	10	0	0	1	10	1	10	0.001	19.6	ACTCCCCTGG	ACTCCCCTGG	11	20	90.9	50.0
};

my ($fh1, $filename1) = File::Temp::tempfile();

print $fh1 $first_file;
close($fh1) || die;

my ($fh2, $filename2) = File::Temp::tempfile();

print $fh2 $second_file;
close($fh2) || die;

my (undef, $filename3) = File::Temp::tempfile(OPEN => 0);

my ($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl');
isnt(Test::Script::Run::last_script_exit_code(), 0, 'Without arguments is should exit with exit code not equals to 1');

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl', ["--query", $filename1] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With file for query is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a subject FASTA file via --subject parameter/, "Test for missing subject file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl', ["--query", $filename1."not_there"] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With non existing file for query is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a query FASTA file via --query parameter/, "Test with non-existing query file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl', ["--subject", $filename2] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With file for subject is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a query FASTA file via --query parameter/, "Test for missing query file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl', ["--query", $filename1, "--subject", $filename2."not_there"] );
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With non existing subject for query is should exit with exit code not equals to 1');
like($stderr, qr/Need to specify a subject FASTA file via --subject parameter/, "Test with non-existing missing subject file");

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl',["--query", $filename1, "--subject", $filename2]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With files for query and subject it should exit with exit code 0');

my $got	= parser($stdout);
my $expected = parser($expected_output);
is_deeply($got,$expected,'output as expected');

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl',["--query", $filename1, "--subject", $filename2, "--out", $filename3]);
is(Test::Script::Run::last_script_exit_code(), 0, 'With files for query, subject, and out it should exit with exit code 0');

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl',["--query", $filename1, "--subject", $filename2, "--out", $filename3]);
isnt(Test::Script::Run::last_script_exit_code(), 0, 'With files for query, subject, and existing out it should exit without exit code 0');

($return,$stdout,$stderr)=run_script('../scripts/MINING_ortholog_blast.pl',["--query", $filename1, "--subject", $filename2, "--out", $filename3, "--overwrite"]);
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
