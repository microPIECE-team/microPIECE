#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin $Script);
use lib "$RealBin/lib/";

use microPIECE;

hello();

sub hello
{
    # version15char is a spacer with 15 character length
    # the current version number is obtained and will be extended to
    # 15 characters in total length to fit into the ASCII art
    my $version15char = $microPIECE::VERSION->normal();
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

__END__

=pod



=cut
