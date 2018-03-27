use strict;
use warnings;

use Test::More;

use FindBin;
use lib "$FindBin::Bin/../scripts/lib";
use Digest::MD5;

BEGIN { use_ok('mining') };

# test fasta import
# cat t/hairpin_high_conf.fa | sed '/^>/s/^\(>[^[:space:]]*\).*$/\1<<</g;' | tr -d "\n" | sed 's/>/\n>/g' | sed '/^$/d;s/<<</\t/;' | sort -k1,1 -k2,2 | sed 's/\t/\n/' | md5sum
#
# 7cf7408ce531fd491e6620a3a8fb1a1e  -

my $dat = mining::parse_fasta("t/hairpin_high_conf.fa");
my $ctx = Digest::MD5->new;
foreach my $header (sort keys %{$dat})
{
    my $fasta_block = ">$header\n".uc($dat->{$header})."\n";
    print STDERR $fasta_block;
    $ctx->add($fasta_block);
}
is($ctx->hexdigest(), "7cf7408ce531fd491e6620a3a8fb1a1e", 'Fasta import can ge exported as expected');

done_testing;
