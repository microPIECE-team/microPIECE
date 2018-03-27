use strict;
use warnings;

use Test::More;

use FindBin;
use lib "$FindBin::Bin/../scripts/lib";
use Digest::MD5;

BEGIN { use_ok('mining') };

# test fasta import
# cat t/hairpin_high_conf.fa | sed '/^>/s/^\(>[^[:space:]]*\).*$/\1<<</g;' | tr -d "\n" | sed 's/>/\n>/g' | sed '/^$/d;s/<<</\t/;' | sort -k2,2 -k1,1 | sed 's/\t/\n/' | md5sum
#
# 0a0a768296513b79e2e1b338044fd134  -

my $dat = mining::parse_fasta("t/hairpin_high_conf.fa");
my $ctx = Digest::MD5->new;
foreach my $header (sort { $dat->{$a} cmp $dat->{$b} || $a cmp $b } keys %{$dat})
{
    my $fasta_block = ">$header\n".uc($dat->{$header})."\n";
    $ctx->add($fasta_block);
}
is($ctx->hexdigest(), "0a0a768296513b79e2e1b338044fd134", 'Fasta import can ge exported as expected');

done_testing;
