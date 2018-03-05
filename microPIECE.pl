#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin $Script);
use lib "$RealBin/lib/";

use IO::Handle;

use microPIECE;
use Log::Log4perl;

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init( \q(
	log4perl.rootLogger                     = INFO, Screen
	log4perl.appender.Screen                = Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr         = 1
	log4perl.appender.Screen.layout         = PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yyyy-MM-dd HH:mm:ss}] %m%n
));

# parse the input parameter
use Getopt::Long;
use Pod::Usage;

# GetOptions
my $opt = {
    version            => undef,
    help               => undef,
    genomeA            => undef,
    genomeB            => undef,
    annotationA        => undef,
    annotationB        => undef,
    clip               => [],
    adapterclip        => undef,
    smallrnaseq        => {},
    adaptersmallrnaseq5=> undef,
    adaptersmallrnaseq3=> undef,
    filterncrnas       => undef,
    config             => undef,
    threads            => 1,
    out                => "out",
    overwrite          => 0,
    scriptdir          => $RealBin."/scripts/",
    testrun            => undef,
    mirna              => undef,
    speciesB_tag       => undef,
};

GetOptions(
    'version|V'            => \$opt->{version},
    'help|h'               => \$opt->{help},
    'genomeA=s'            => \$opt->{genomeA},
    'genomeB=s'            => \$opt->{genomeB},
    'gffA|annotationA=s'   => \$opt->{annotationA},
    'gffB|annotationB=s'   => \$opt->{annotationB},
    'clip=s@'              => \$opt->{clip},
    'adapterclip=s'        => \$opt->{adapterclip},
    'smallrnaseq=s%'       => sub { push(@{$opt->{smallrnaseq}{$_[1]}}, $_[2]) },
    'adaptersmallrnaseq5=s'=> \$opt->{adaptersmallrnaseq5},
    'adaptersmallrnaseq3=s'=> \$opt->{adaptersmallrnaseq3},
    'filterncrnas=s'       => \$opt->{filterncrnas},
    'config=s'             => \$opt->{config},
    'threads=i'            => \$opt->{threads},
    'overwrite'            => \$opt->{overwrite},
    'testrun'              => \$opt->{testrun},
    'out=s'                => \$opt->{out},
    'mirna=s'              => \$opt->{mirna},
    'speciesBtag=s'           => \$opt->{speciesB_tag},
    ) || pod2usage(1);

# split clip files if required
$opt->{clip} = [ split(",", join(",", @{$opt->{clip}})) ];
# split rnaseq files if required
foreach my $cond (keys %{$opt->{smallrnaseq}})
{
    $opt->{smallrnaseq}{$cond} = [ split(",", join(",", @{$opt->{smallrnaseq}{$cond}})) ];
}

# help
$opt->{help} && pod2usage(1);

# version
if($opt->{version}){
	print $microPIECE::VERSION->normal(),"\n";
	exit 0;
}

microPIECE::hello();
STDOUT->flush();

microPIECE::check_requirements($opt);

microPIECE::print_settings($opt);

microPIECE::run_mining($opt);

microPIECE::run_clip($opt);

microPIECE::run_targetprediction($opt);

microPIECE::transfer_resultfiles($opt);

__END__

=pod

=head1 NAME

microPIECE - microRNA pipeline enhanced by CLIP experiments

=head1 SYNOPSIS

./microPIECE.pl \
  --genomeA <speciesA_genome.fa> \
  --genomeB <speciesB_genome.fa> \
  --annotationA <speciesA_genome.gff> \
  --annotationB <speciesB_genome.gff> \
  --clip <ago_clip_seq.fastq> \
  --mirnas <mature_miRNAs.fa> \

=head1 DESCRIPTION

The microPIECE (microRNA pipeline enhanced by CLIP experiments) 
takes the AGO-CLIP data from a speciesA and transfers it to a speciesB. 
Given a set of miRNAs from speciesB it then predicts their targets on the transfered CLIP regions.

For the minimal workflow it needs a genome file, as well as its annotation file in GFF format for speciesA and speciesB. 
For speciesA at least one AGO-CLIP dataset is needed and speciesB needs a set of miRNAs for the target prediction.
For the full workflow, a set of smallRNA-sequencing data is additionally needed and a set of non-coding RNAs can be 
provided as filter. The pipeline uses the smallRNA data for the mining of novel microRNAs and the completion of 
the given miRNA dataset, if needed. It further performs expression calculation, isoform detection, genomic loci 
identification and orthology determination.

=head1 EXAMPLE

./microPIECE.pl   \
  --genomeA testset/NC_035109.1_reduced_AAE_genome.fa  \
  --genomeB testset/NC_007416.3_reduced_TCA_genome.fa   \
  --annotationA testset/NC_035109.1_reduced_AAE_genome.gff   \
  --annotationB testset/NC_007416.3_reduced_TCA_genome.gff   \
  --clip testset/SRR5163632_aae_clip_reduced.fastq,testset/SRR5163633_aae_clip_reduced.fastq,testset/SRR5163634_aae_clip_reduced.fastq   \
  --clip testset/SRR5163635_aae_clip_reduced.fastq,testset/SRR5163636_aae_clip_reduced.fastq,testset/SRR5163637_aae_clip_reduced.fastq --adapterclip GTGTCAGTCACTTCCAGCGG  \
  --overwrite \
  --smallrnaseq a=testset/tca_smallRNAseq_rna_contaminated.fastq \
  --adaptersmallrnaseq3=TGGAATTCTCGGGTGCCAAGG \
  --adaptersmallrnaseq5 GTTCAGAGTTCTACAGTCCGACGATC \
  --filterncrnas testset/TCA_all_ncRNA_but_miR.fa \
  --speciesB tca 2>&1 | tee out.log

=head1 OUTPUT

=head2 Output data

=over 4

=item mature miRNA set: C<mature_combined_mirbase_novel.fa>

mature microRNA set, containing novels and miRBase-completed (if mined), together with the known miRNAs from miRBase

=item precursor miRNA set: C<hairpin_combined_mirbase_novel.fa>

precursor microRNA set, containing novels (if mined), together with the known miRNAs from miRBase

=item mature miRNA expression per condition: C<miRNA_expression.csv>

Semicolon-separated file --> rpm;condition;miRNA

=item miRDeep2 mining result in HTML/CSV C<mirdeep_output.html/csv>

the standard output HTML/CSV file of miRDeep2

=item all library support-level target predictions: C<*_miranda_output.txt>

miranda output, reduced to the lines, starting with > only

=item all library support-level CLIP transfer .bed files: C<*transfered_merged.bed>

bed-file of the transferred CLIP-regions in speciesB transcriptome

=back

=head1 HISTORY
1.0.0 - March, 2nd 2018 : Initial release - no known bugs.

=head2 INPUT PARAMETERS:

=over 4

=item --version|-v

version of this pipeline

=item --help|-h

prints a helpful help message

=item --genomeA

Genome of the species with the CLIP data

=item --genomeB

Genome of the species where we want to predict the miRNA targets

=item --gffA

GFF annotation of speciesA

=item --gffB

GFF annotation of speciesB

=item --clip

Comma-separated CLIP-seq .fastq files in Format

		--clip con1_rep1_clip.fq,con1_rep2_clip.fq,con2_clip.fq
		OR
		--clip con1_rep1_clip.fq --clip con1_rep2_clip.fq --clip con2_clip.fq

=item 	--adapterclip

Sequencing-adapter of CLIP reads

=item --smallrnaseq

Comma-separated smallRNA-seq .fastq files, initialized with 'condition=' in Format

		--smallrnaseq con1=A.fastq,B.fastq --smallrnaseq con2=C.fq
		OR
		--smallrnaseq con1=A.fastq --smallrnaseq con1=B.fastq --smallrnaseq con2=C.fq

=item --adaptersmallrnaseq5

5' adapter of smallRNA-seq reads

=item --adaptersmallrnaseq3

3' adapter of smallRNA-seq reads

=item --filterncrnas

Multi-fasta file of ncRNAs to filter smallRNA-seq reads

=item --threads

Number of threads to be used

=item --overwrite

set this parameter to overwrite existing files

=item --testrun

sets this pipeline to testmode (accounting for small testset in piranha)

=item --out

output folder

=item --mirnas

miRNA set, if set, mining is disabled and this set is used for prediction

=item --speciesBtag

Three letter code of speciesB

=back

=cut
