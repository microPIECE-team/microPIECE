#! /usr/bin/perl
# Author : Daniel Amsel - daniel.amsel@ime.fraunhofer.de
use strict;
use warnings;
use Getopt::Long;
use IPC::Run3;


my $query;
my $subject;
my $out = "-";
my $overwrite;
my $threads=1;

my $allowed_gap_mismatches = 1;  # number of mismatches/gaps allowed in alignment (if not in first 10 positions)

GetOptions(
    'query=s' => \$query,
    'subject=s' => \$subject,
    'out=s' => \$out,
    'force|overwrite' => \$overwrite,
    'threads=i' => \$threads
    );

die "Need to specify a query FASTA file via --query parameter\n" unless (defined $query && -e $query);
die "Need to specify a subject FASTA file via --subject parameter\n" unless (defined $subject && -e $subject);
die "Output file exists and --force/--overwrite was not specified\n" if ($out ne "-" && (-e $out) && (! $overwrite));

# generate database
my $cmd = ["makeblastdb", "-in", $subject, "-dbtype", "nucl"];
run3($cmd, undef, \(my $stdout), \(my $stderr));
if ($? != 0)
{
    die("Error calling command: ".join(" ", @{$cmd})."\n***************\nSTDERR\n***************\n$stderr\n\n***************\nSTDOUT\n***************\n$stdout\n");
}

my $blastcmd = "blastn -db $subject -query $query -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen slen' -num_threads $threads -word_size 4 -evalue 10000 -strand plus";
open(BLAST, $blastcmd."|") || die "Unable to open pipe via command '$blastcmd': $!\n";

my $fh;
if ($out ne "-")
{
    open($fh, ">", $out) || die "Unable to open file '$out': $!\n";
} else {
    $fh = *STDOUT;
}

while (<BLAST>)
{
    chomp;
    my %fields = ();
    @fields{qw(query subject identity length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen slen)} = split("\t", $_);

    # need to have minimum alignment length of 10
    next if ($fields{length} < 10);

    # alignments should start on query and on subject from first
    # nucleotid
    next unless (($fields{qstart}==1) && ($fields{sstart}==1));

    # the seed region (first 10 bp) have to be aligned without
    # mismatch or gap
    my @qaln = split("",$fields{qseq});
    my @saln = split("",$fields{sseq});

    my $bool_keep = 1;	# assume to keep the BLAST hit
    for (my $i = 0; $i<10; $i++)
    {
	if($qaln[$i] ne $saln[$i])
	{
	    $bool_keep = 0;
	    last; # no need to run through all 10 bp if we found a
		  # mismatch/gap
	}
    }

    next unless ($bool_keep);

    # last criterion for acceptance is: rest of alignment must not
    # contain more than one mismatch/gap
    my $num_gap_mismatches = 0;
    for(my $i = 10; $i < int(@qaln); $i++)
    {
	# check for gap or mismatch
	if (($qaln[$i] eq "-" || $saln[$i] eq "-") || ($qaln[$i] ne $saln[$i]))
	{
	    $num_gap_mismatches++;
	}

	if ($num_gap_mismatches > $allowed_gap_mismatches)
	{
	    $bool_keep = 0;
	    last;
	}
    }

    next unless ($bool_keep);

    # calulate the coverage for query and subject
    # we need the aligned sequence length for query/subject without
    # gaps
    my $qseq_wo_gaps = $fields{qseq};
    $qseq_wo_gaps =~ tr/-//;

    my $sseq_wo_gaps = $fields{sseq};
    $sseq_wo_gaps =~ tr/-//;

    my $qcoverage = sprintf("%.1f", (length($qseq_wo_gaps)/$fields{qlen}*100));
    my $scoverage = sprintf("%.1f", (length($sseq_wo_gaps)/$fields{slen}*100));

    print $fh join("\t", ($_, $qcoverage, $scoverage)), "\n";

}
if ($out eq "-")
{
    close($fh) || die "Unable to close file '$out': $!\n";
}
close(BLAST) || die "Unable to close pipe for command '$cmd': $!\n";
