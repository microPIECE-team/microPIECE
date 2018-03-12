package bed;

use strict;
use warnings;

sub parser{
    my $p_stdout		= $_[0];

    my @p_stdout_array	= ();

    foreach my $bed_line (split("\n",$p_stdout)){
	next if ($bed_line =~ /^#/);

	my @bed_fields		= split("\t",$bed_line);
	my @key_value_pairs 	= split(";",$bed_fields[3]);

	$bed_fields[3]		= {};

	foreach(@key_value_pairs){
	    my ($key,$value)= split("=",$_,2);
	    $bed_fields[3]{$key}=$value;
	}

	push(@p_stdout_array,\@bed_fields);
    }

    return( [ sort
	      {
		  $a->[0] cmp $b->[0] ||
		      $a->[1] <=> $b->[1] ||
		      $a->[2] <=> $b->[2] ||
		      $a->[5] cmp $b->[5] ||
		      $a->[3]{counts} cmp $b->[3]{counts}
	      } @p_stdout_array ] );
}


1;
