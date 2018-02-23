#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin $Script);
use lib "$RealBin/lib/";

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
    adaptersmallrnaseq => undef,
    filterncrnas       => undef,
    config             => undef,
    threads            => 1,
    out                => "out",
    overwrite          => 0,
    scriptdir          => $RealBin."/scripts/",
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
    'adaptersmallrnaseq=s' => \$opt->{adaptersmallrnaseq},
    'filterncrnas=s'       => \$opt->{filterncrnas},
    'config=s'             => \$opt->{config},
    'threads=i'            => \$opt->{threads},
    'overwrite'            => \$opt->{overwrite},
    'out=s'                => \$opt->{out}
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

microPIECE::print_settings($opt);

microPIECE::check_requirements($opt);

microPIECE::print_settings($opt);

microPIECE::run_mining($opt);

microPIECE::run_clip($opt);

microPIECE::run_targetprediction($opt);

__END__

=pod



=cut
