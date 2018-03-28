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

$dat = mining::parse_mirbase_dat("t/miRNA_high_conf.dat", "");

@{$dat} = sort {$a->{seq} cmp $b->{seq} || $a->{precursor} cmp $b->{precursor} } @{$dat};
$ctx = Digest::MD5->new;
foreach my $entry (@{$dat})
{
    my $fasta_block = ">".$entry->{precursor}."\n".uc($entry->{seq})."\n";
    $ctx->add($fasta_block);
}
is($ctx->hexdigest(), "0a0a768296513b79e2e1b338044fd134", 'mirbase *.dat import contains all haipins');

# mature checksum:
# cat t/mature_high_conf.fa | sed '/^>/s/^\(>[^[:space:]]*\).*$/\1<<</g;' | tr -d "\n" | sed 's/>/\n>/g' | sed '/^$/d;s/<<</\t/;' | sort -k2,2 -k1,1 | sed 's/\t/\n/' | md5sum
#
# fb0a6f1793a8c284df4834c6cfbbe3e4  -
my @tmp = ();
foreach my $entry (@{$dat})
{
    foreach my $mature (@{$entry->{matures}})
    {
	push(@tmp, { name => $mature->{name}, seq => substr($entry->{seq}, $mature->{start}-1, $mature->{stop}-$mature->{start}+1) });
    }
}
@tmp = sort { $a->{seq} cmp $b->{seq} || $a->{name} cmp $b->{name} } (@tmp);

$ctx = Digest::MD5->new;
my %seen = ();
foreach my $entry (@tmp)
{
    my $fasta_block = ">".$entry->{name}."\n".uc($entry->{seq})."\n";
    $ctx->add($fasta_block) unless (exists $seen{$fasta_block});
    $seen{$fasta_block}++;
}
is($ctx->hexdigest(), "fb0a6f1793a8c284df4834c6cfbbe3e4", 'mirbase *.dat import contains all matures');

@tmp = ();
foreach my $entry (@{$dat})
{
    foreach my $mature (@{$entry->{matures}})
    {
	push(@tmp, { name => $mature->{name}, seq => $mature->{seq} });
    }
}
@tmp = sort { $a->{seq} cmp $b->{seq} || $a->{name} cmp $b->{name} } (@tmp);

$ctx = Digest::MD5->new;
%seen = ();
foreach my $entry (@tmp)
{
    my $fasta_block = ">".$entry->{name}."\n".uc($entry->{seq})."\n";
    $ctx->add($fasta_block) unless (exists $seen{$fasta_block});
    $seen{$fasta_block}++;
}
is($ctx->hexdigest(), "fb0a6f1793a8c284df4834c6cfbbe3e4", 'mirbase *.dat import contains all matures');

done_testing;
