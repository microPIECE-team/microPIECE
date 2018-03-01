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
    speciesBtag        => undef,
};
# INPUT PARAMETERS:
#	--version|-v		:= version of this pipeline
#	--help|-h		:= prints a helpful help message
#	--genomeA		:= Genome of the species with the CLIP data
#	--genomeB		:= Genome of the species where we want to predict the miRNA targets
#	--gffA			:= GFF annotation of speciesA
#	--gffB			:= GFF annotation of speciesB
#	--clip			:= Comma-separated CLIP-seq .fastq files
	#	--clip con1_rep1_clip.fq,con1_rep2_clip.fq,con2_clip.fq
	#	OR
	#	--clip ron1_rep1_clip.fq --clip con1_rep2_clip.fq --clip con2_clip.fq
#	--adapterclip		:= Sequencing-adapter of CLIP reads
#	--smallrnaseq		:= Comma-separated smallRNA-seq .fastq files, initialized with 'condition='
	#	--smallrnaseq con1=A.fastq,B.fastq --smallrnaseq con2=C.fq
	#	OR
	#	--smallrnaseq con1=A.fastq --smallrnaseq con1=B.fastq --smallrnaseq con2=C.fq
#	--adaptersmallrnaseq5	:= 5' adapter of smallRNA-seq reads
#	--adaptersmallrnaseq3	:= 3' adapter of smallRNA-seq reads
#	--filterncrnas		:= Multi-fasta file of ncRNAs to filter smallRNA-seq reads
#	--threads		:= Number of threads to be used
#	--overwrite		:= set this parameter to overwrite existing files 
#	--testrun		:= sets this pipeline to testmode (accounting for small testset in piranha)
#	--out			:= output folder
#	--mirnas		:= miRNA set, if set, mining is disabled and this set is used for prediction
#	--speciesBtag		:= 3letter code of speciesB
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
    'speciesBtag=s'        => \$opt->{speciesBtag},
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

__END__

=pod



=cut
