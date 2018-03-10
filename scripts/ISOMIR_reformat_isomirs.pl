#! /usr/bin/perl
use strict;
use warnings;


my $input_list		= $ARGV[0];	# comma separated replicates
my $condition_ID	= $ARGV[1];	# from db

my @input_array		= split(",",$input_list);
my %condition_hash;	#{condition_file}	= \%rep_hash;

my @required_fields = qw(seq freq mir mism add t5 t3 ambiguity);

my %rep_hash;	 #{miRBaseID;mism;add;t5;t3;ambiguity} = freq/ambiguous

foreach my $input_file	(@input_array){
    open(TMP,"<",$input_file) || die "Unable to open input file '$input_file': $!\n";
    
    my $fieldlist	= <TMP>;
    chomp($fieldlist);
    my @possible_fields = split(/\s+/, $fieldlist);

    # check if required fields are existing
    foreach my $req_field (@required_fields)
    {
	unless (grep {$_ eq $req_field} (@possible_fields))
	{
	    die("Missing field '$req_field' in input file '$input_file'\n");
	}
    }
    
    my $tmp_total_read_count= 0;
    while(<TMP>)
    {
	chomp;
	my %fields = ();
	@fields{@possible_fields} = split("\t", $_);

	# check if one of the mism add t5 t3 fields is 0
	foreach my $field (qw(mism add t5 t3))
	{
	    if ($fields{$field} eq "0")
	    {
		$fields{$field} = 'NULL';
	    }
	}

	my $div_c = $fields{freq}/$fields{ambiguity};
	
	$tmp_total_read_count += $div_c;	# count all reads, divided by ambiguity, but also the non-isoforms
	
	next if(($fields{mism} eq 'NULL') && ($fields{add} eq 'NULL') && ($fields{t5} eq 'NULL') && ($fields{t3} eq 'NULL'));

	# generate a key from the fields: mir mism add t5 t3 seq
	my $key = join(";", map {$fields{$_}} (qw(mir mism add t5 t3 seq)));

	$rep_hash{$key}{fields} = \%fields;
	$rep_hash{$key}{rpm}{$input_file} = $div_c;
    }

    close(TMP) || die "Unable to close input file '$input_file': $!\n";

    # calc RPM
    my %rpm_hash;
    foreach my $key (keys %rep_hash){
	next unless (exists $rep_hash{$key}{rpm}{$input_file}); 

	$rep_hash{$key}{rpm}{$input_file} = $rep_hash{$key}{rpm}{$input_file}/$tmp_total_read_count*1000000;
    }
}

my $condition_number = scalar(@input_array); # number of conditions : e.g. 4

# calculate the mean for each rpm
foreach my $key (sort keys %rep_hash)
{
    my $mean = 0;
    foreach my $condition (keys %{$rep_hash{$key}{rpm}})
    {
	$mean += $rep_hash{$key}{rpm}{$condition};
    }
    $mean = $mean/$condition_number;

    print join(";", (
		   (map {$rep_hash{$key}{fields}{$_}} (qw(mir mism add t5 t3 seq))),
		   $mean,
		   $condition_ID
	       )
	), "\n";
}
