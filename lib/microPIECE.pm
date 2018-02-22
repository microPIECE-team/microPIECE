package microPIECE;

use version 0.77; our $VERSION = version->declare("v0.9.0");

use Log::Log4perl;
use Data::Dumper;

=pod

Just a little welcome screen indicating our current version number

=cut

sub hello
{
    # version15char is a spacer with 15 character length
    # the current version number is obtained and will be extended to
    # 15 characters in total length to fit into the ASCII art
    my $version15char = $VERSION->normal();
    my $required_spaces = 15-length($version15char)-1;
    $version15char .= " "x$required_spaces;
    print "
                    @ @ @ @ @ @
                    @       @ @
                    @ @   @ @ @
          @ @ @ @ @ @ @   @ @ @ @ @ @ @
          @ @                         @
          @ @                         @
  @ @ @   @ @                     @ @ @
@ @   @ @ @ @    m i c r o      @ @
@       @ @ @                   @
@                P I E C E      @
@       @ @ @                   @ @
@ @ @ @ @ @ @    $version15char   @ @ @
    @     @ @                         @
          @ @                         @
          @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
                  @ @         @ @
                  @ @         @ @
                  @ @ @ @ @ @ @ @
                  @ @         @ @
                  @ @         @ @
                  @ @ @ @ @ @ @ @
                  @ @         @ @
                  @ @         @ @
                  @ @ @ @ @ @ @ @
                  @ @         @ @
                  @ @         @ @

"
}

=pod

Checking requirements for the pipeline

=cut

sub check_requirements {
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    # first we check if mandatory parameters have been specified
    unless (exists $opt->{genomeA} && -e $opt->{genomeA})
    {
	$L->logdie("Missing parameter for --genomeA or file is inaccessable\n");
    }

    unless (exists $opt->{genomeB} && -e $opt->{genomeB})
    {
	$L->logdie("Missing parameter for --genomeB or file is inaccessable\n");
    }
}

=pod

Print settings

=cut

sub print_settings {

    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    $L->info(Dumper($opt));
}

=pod

Running the mining part of the microPIECE

=cut

sub run_mining {

    my ($opt) = @_;

}

=pod

Running the CLIP part of the microPIECE

=cut

sub run_clip {

    my ($opt) = @_;

}

=pod

Running the target prediction part of the microPIECE

=cut

sub run_targetprediction {

    my ($opt) = @_;

}


1;
