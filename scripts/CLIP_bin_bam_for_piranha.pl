#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my @bamfiles = ();
my $binsize  = 20;

GetOptions(
    'b|bam=s@' => \@bamfiles,
    's|size=i' => \$binsize,
    ) || die;

# definde the bams if comma seperated
@bamfiles = split(",", join(",", @bamfiles));

# try to call samtools
my $samtools = qx(samtools --version);
if ($? != 0)
{
    die "Wrong version or no samtools available\n";
}

foreach my $bam (@bamfiles)
{
    my $cmd = "samtools view -H $bam";
    my $header = qx($cmd);
    if ($? != 0)
    {
	die "Error running command '$cmd'\n";
    }
}
