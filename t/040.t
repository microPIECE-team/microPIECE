#!/usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

my($return,$stdout,$stderr)=run_script('../scripts/046_merge_bed_files.pl');
like($stderr,qr/Example defining to classes with two files each and output to merged.bed/,'print help when no args provided');

($return,$stdout,$stderr)=run_script('../scripts/046_merge_bed_files.pl',[qw(--input bed1=040/041_testset_01.bed)]);
#is($Test::Script::Run::last_script_exit_code,0,'script with one bed file as input');
#print "$return\n$stdout\n$stderr\n";


my $got	= &parser($stdout);
my $expected = &parser(join('',<DATA>));	
is_deeply($got,$expected,'output as expected');

done_testing();


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
	return(\@p_stdout_array);
}
__DATA__
NC_010241.1	0	1	length=1;counts=1/1	.	+
NC_010241.1	30	45	length=15;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,2/2,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	10	20	length=10;counts=2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2,2/2	.	+
NC_010241.1	10	20	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	-
NC_010241.1	50	60	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	60	65	lenght=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+
NC_010241.1	70	80	length=10;counts=1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1,1/1	.	+
