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
	log4perl.appender.Screen.layout.ConversionPattern = [%d{yy-MM-dd HH:mm:ss}] %m%n
));

# parse the input parameter
use Getopt::Long;
use Pod::Usage;

# GetOptions
my %opt = ();

GetOptions( # use %opt as defaults (demo-set)
	    \%opt, qw(
              version|V
              help|h
	)
    ) || pod2usage(1);

# help
$opt{help} && pod2usage(1);

# version
if($opt{version}){
	print $microPIECE::VERSION->normal(),"\n";
	exit 0;
}

microPIECE::hello();

microPIECE::check_dependencies();

microPIECE::run_mining(\%opt);

__END__

=pod



=cut
