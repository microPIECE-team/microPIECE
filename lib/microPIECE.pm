package microPIECE;

use strict;
use warnings;

use version 0.77; our $VERSION = version->declare("v0.9.0");

use Log::Log4perl;
use Data::Dumper;
use Cwd;
use File::Path;
use File::Basename;

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

    $opt->{basedir} = getcwd."/";

    # first we check if mandatory parameters have been specified
    unless (exists $opt->{genomeA} && -e $opt->{genomeA})
    {
	$L->logdie("Missing parameter for --genomeA or file is inaccessable\n");
    }
    $opt->{genomeA} = $opt->{basedir}.$opt->{genomeA};

    unless (exists $opt->{genomeB} && -e $opt->{genomeB})
    {
	$L->logdie("Missing parameter for --genomeB or file is inaccessable\n");
    }
    $opt->{genomeB} = $opt->{basedir}.$opt->{genomeB};

    unless (exists $opt->{annotationA} && -e $opt->{annotationA})
    {
	$L->logdie("Missing parameter for --annotationA or file is inaccessable\n");
    }
    $opt->{annotationA} = $opt->{basedir}.$opt->{annotationA};

    unless (exists $opt->{annotationB} && -e $opt->{annotationB})
    {
	$L->logdie("Missing parameter for --annotationB or file is inaccessable\n");
    }
    $opt->{annotationB} = $opt->{basedir}.$opt->{annotationB};

    if (exists $opt->{clip} && @{$opt->{clip}} > 0)
    {
	for( my $i=0; $i<@{$opt->{clip}}; $i++ )
	{
	    $L->logdie("Missing parameter for --clip or file is inaccessable\n") unless (-e $opt->{clip}[$i]);
	    $opt->{clip}[$i] = $opt->{basedir}.$opt->{clip}[$i];
	}
    } else {
	$L->logdie("Missing parameter for --clip or file is inaccessable\n")
    }

    unless (exists $opt->{adapterclip} && length($opt->{adapterclip})>0 && $opt->{adapterclip} =~ /^[ACGT]+$/i)
    {
	$L->logdie("Missing parameter for --adapterclip or unexpected characters provided");
    }

    if (-e $opt->{out})
    {
	unless ( $opt->{overwrite} )
	{
	    $L->error(sprintf("Output folder %s exists. Use --overwrite to overwrite it", $opt->{out}));
	    exit;
	} else {
	    rmtree($opt->{out}, {error => \my $err} );
	    if (@$err) {
		for my $diag (@$err) {
		    my ($file, $message) = %$diag;
		    if ($file eq '') {
			$L->logdie("general error: $message");
		    }
		    else {
			$L->logdie("problem unlinking $file: $message");
		    }
		}
	    }
	}
    }

    mkdir($opt->{out}) || $L->logdie("Unable to create output file: $!");

    $opt->{basedir} = $opt->{basedir}.$opt->{out}."/";
    chdir($opt->{basedir});

    # we need to run clip
    $opt->{run_clip} = 1;
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

    my $L = Log::Log4perl::get_logger();

    unless (exists $opt->{run_mining} && $opt->{run_mining})
    {
	$L->info("Skipping mining step due to missing parameters");
	return;
    }

    $L->info("Starting mining step");

    $L->info("Finished mining step");

}

=pod

Running the CLIP part of the microPIECE

=cut

sub run_clip {

    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    unless (exists $opt->{run_clip} && $opt->{run_clip})
    {
	$L->info("Skipping CLIP step due to missing parameters");
	return;
    }

    $L->info("Starting CLIP step");

    run_proteinortho($opt);
    run_CLIP_adapter_trimming($opt);
    run_CLIP_build_db($opt);
    run_CLIP_mapping($opt);
    run_CLIP_piranha($opt);
    run_CLIP_bedtools_merge($opt);
    run_CLIP_filterbed($opt);
    run_CLIP_clip_mapper($opt);
    run_CLIP_process($opt);
    run_CLIP_transfer($opt);

    $L->info("Finished CLIP step");

}

sub run_CLIP_transfer
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my @inputfiles = glob("clip_merged_*of*BEDfilter_mapGFF_minLen*_min*_max*_sort_UC.fasta");

    my $file_uniqueA     = $opt->{basedir}.basename($opt->{annotationA}, ".gff")."_unique.csv";
    my $file_uniqueA_log = $opt->{basedir}.basename($opt->{annotationA}, ".gff")."_unique.err";
    my @cmd = ($opt->{scriptdir}."085_parse_gff_return_longest_transcript.pl", $opt->{annotationA}, ">", $file_uniqueA, "2>", $file_uniqueA_log);
    run_cmd($L, \@cmd);

    my $file_uniqueB     = $opt->{basedir}.basename($opt->{annotationB}, ".gff")."_unique.csv";
    my $file_uniqueB_log = $opt->{basedir}.basename($opt->{annotationB}, ".gff")."_unique.err";
    @cmd = ($opt->{scriptdir}."085_parse_gff_return_longest_transcript.pl", $opt->{annotationB}, ">", $file_uniqueB, "2>", $file_uniqueB_log);
    run_cmd($L, \@cmd);

    foreach my $file (@inputfiles)
    {
	my $needle_csv = $opt->{basedir}.basename($file, ".fasta")."_needle.csv";
	my $needle_aln = $opt->{basedir}.basename($file, ".fasta")."_needle.aln";
	@cmd = ($opt->{scriptdir}."081_map_clip_gff_needle.pl", $file_uniqueA, $opt->{proteinortho}, $file, $file_uniqueB, $opt->{genomeB}, $needle_csv, ">", $needle_aln);
	run_cmd($L, \@cmd);
    }
}

sub run_CLIP_process
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my $min = 22;
    my $max = 50;

    my @inputfiles = glob("clip_merged_*of*BEDfilter_mapGFF_minLen0.bed");

    foreach my $file (@inputfiles)
    {
	my $sorted_bed = sprintf("%s%s_min%i_max%i_sort.bed",      $opt->{basedir}, basename($file, ".bed"), $min, $max);
	my $fasta =      sprintf("%s%s_min%i_max%i_sort.fasta",    $opt->{basedir}, basename($file, ".bed"), $min, $max);
	my $fastaUC =    sprintf("%s%s_min%i_max%i_sort_UC.fasta", $opt->{basedir}, basename($file, ".bed"), $min, $max);
	my @cmd = ($opt->{scriptdir}."071_bedtool_discard_sizes.pl", $file, $min, $max);
	my $output = run_cmd($L, \@cmd);

	# sort the file
	my @dat = ();
	foreach my $line (split(/\n/, $output))
	{
	    my @fields = split("\t", $line);
	    push(@dat, \@fields);
	}

	@dat = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[2] } (@dat);

	open(FH, ">", $sorted_bed) || $L->logdie("Unable to open '$sorted_bed': $!");
	foreach my $fields (@dat)
	{
	    print FH join("\t", @{$fields}), "\n";
	}
	close(FH) || $L->logdie("Unable to close '$sorted_bed': $!");

	# run bedtools getfasta
	@cmd = ("bedtools", "getfasta", "-s", "-name", "-fi", $opt->{genomeA}, "-bed", $sorted_bed, );
	my $fastaoutput = run_cmd($L, \@cmd);

	open(FH, ">", $fasta) || $L->logdie("Unable to open '$fasta': $!");
	foreach my $line (split(/\n/, $fastaoutput))
	{
	    # bedtools
	    if ($line =~ /^>/ && $line =~ s/\([-+]\)$//)
	    {
		$L->warn("Found (+/-) as last information in fasta header. Assuming addition through 'bedtools getfasta' and removing that");
	    }

	    print FH $line, "\n";
	}
	close(FH) || $L->logdie("Unable to close '$fasta': $!");

	# convert fasta to upper case
	@cmd = ($opt->{scriptdir}."072_fasta_uc_and_filter4annotations.pl", $fasta, ">", $fastaUC);
	run_cmd($L, \@cmd);
    }
}

sub run_CLIP_clip_mapper
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my $minlength = 0;

    my @inputfiles = glob("clip_merged_*of*BEDfilter.bed");

    foreach my $file (@inputfiles)
    {
	my $outputname = $opt->{basedir}.basename($file, ".bed")."_mapGFF_minLen0.bed";
	my @cmd = ($opt->{scriptdir}."051_clip_mapper.pl", $file, $opt->{annotationA}, $minlength, ">", $outputname);
	run_cmd($L, \@cmd);
    }
}

sub run_CLIP_filterbed
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my @cmd = ($opt->{scriptdir}."049_bed2signal.pl", $opt->{basedir}."clip_merged.bed");

    my $num_fields = int(@{$opt->{clip}});

    for(my $i=1;$i<=@{$opt->{clip}};$i++)
    {
	my $filtered_bed_out = sprintf("%sclip_merged_%dof%dBEDfilter.bed", $opt->{basedir}, $i, $num_fields);
	run_cmd($L, [@cmd, $i, ">", $filtered_bed_out] );
    }
}


sub run_CLIP_bedtools_merge
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my @cmd = ($opt->{scriptdir}."046_merge_bed_files.pl", "--output", $opt->{basedir}."clip_merged.bed", "--log", "merging_bed_files.log");

    for(my $i=0;$i<@{$opt->{clip}};$i++)
    {
	my $sortedpiranhafile = $opt->{basedir}.basename($opt->{clip}[$i]).".piranha.sorted.bed";
	push(@cmd, ("--input", basename($opt->{clip}[$i])."=".$sortedpiranhafile));
    }

    run_cmd($L, \@cmd);
}

sub run_CLIP_piranha
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    foreach my $clipfile (@{$opt->{clip}})
    {
	my $bedfile           = $opt->{basedir}.basename($clipfile).".bed";
	my $piranhafile       = $opt->{basedir}.basename($clipfile).".piranha.bed";
	my $sortedpiranhafile = $opt->{basedir}.basename($clipfile).".piranha.sorted.bed";
	my @cmd = ("Piranha", "-o", $piranhafile, "-s", $bedfile);
	if (exists $opt->{testrun} && $opt->{testrun})
	{
	    $L->info("TESTRUN was activated though --testrun option. This increases the p-value threshold for Piranha to 20%!!! Please use only for the provided testset and NOT(!!!) for real analysis!!!");
	    push(@cmd, ("-p", 0.2));
	}
	run_cmd($L, \@cmd);

	# own sort routine, was originally based on a sort call,
	# nevertheless, I want to scan for lines containing -nan from
	# Piranha output and skip that lines, due to they will produce
	# errors later.
	my @dat = ();
	open(FH, "<", $piranhafile) || $L->logdie("Unable to open '$piranhafile': $!");
	while (<FH>)
	{
	    chomp;
	    my @fields = split("\t", $_);
	    if ($fields[-1] eq "-nan")
	    {
		$L->warn(sprintf("Line %d from file '%s' contained '-nan' and will be skipped", $., $piranhafile));
	    } else {
		push(@dat, \@fields);
	    }
	}
	close(FH) || $L->logdie("Unable to close '$piranhafile': $!");
	@dat = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[2] } (@dat);

	open(FH, ">", $sortedpiranhafile) || $L->logdie("Unable to open '$sortedpiranhafile': $!");
	foreach my $fields (@dat)
	{
	    print FH join("\t", @{$fields}), "\n";
	}
	close(FH) || $L->logdie("Unable to close '$sortedpiranhafile': $!");
    }
}

sub run_CLIP_mapping
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    foreach my $clipfile (@{$opt->{clip}})
    {
	my $trimmedfile = $opt->{basedir}.basename($clipfile).".trim";
	my $bamfile     = $opt->{basedir}.basename($clipfile).".bam";
	my $bedfile     = $opt->{basedir}.basename($clipfile).".bed";
	my @cmd = ("gsnap", "-N", 1, "-B", 5, "-O", "-A", "sam", "-t", $opt->{threads}, "-D", "speciesA_db", "-d", "speciesA", $trimmedfile, "|", "samtools", "view", "-Sb", "-", "|", "samtools", "sort", "-o", $bamfile, "-");
	run_cmd($L, \@cmd);
	@cmd = ("bedtools", "bamtobed", "-i", $bamfile, ">", $bedfile);
	run_cmd($L, \@cmd);
    }
}

sub run_CLIP_build_db
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();
    my @cmd = ("gmap_build", "-D", "speciesA_db", "-k", 15, "-d", "speciesA", $opt->{genomeA});
    run_cmd($L, \@cmd);
}


sub run_CLIP_adapter_trimming
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    foreach my $clipfile (@{$opt->{clip}})
    {
	my $outfile = $opt->{basedir}.basename($clipfile).".trim";
	my @cmd = ("cutadapt", "-a", $opt->{adapterclip}, "-m", 20, "--trim-n", "-o", $outfile, $clipfile);
	run_cmd($L, \@cmd);
    }
}

=pod

Run proteinortho for ortholog detection

=cut

sub run_proteinortho
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    # extract proteins from annotation and genome file
    $opt->{proteinA} = $opt->{basedir}."/"."proteinA.fa";
    $opt->{proteinB} = $opt->{basedir}."/"."proteinB.fa";
    my @cmd = ("gffread", $opt->{annotationA}, "-y", $opt->{proteinA}, "-F", "-g", $opt->{genomeA});
    run_cmd($L, \@cmd);
    @cmd = ("gffread", $opt->{annotationB}, "-y", $opt->{proteinB}, "-F", "-g", $opt->{genomeB});
    run_cmd($L, \@cmd);

    # create blast databases
    @cmd = ("makeblastdb", "-in", $opt->{proteinA}, "-dbtype", "prot");
    run_cmd($L, \@cmd);
    @cmd = ("makeblastdb", "-in", $opt->{proteinB}, "-dbtype", "prot");
    run_cmd($L, \@cmd);

    # run proteinortho
    @cmd = ("proteinortho5.pl", "-clean", "-project=microPIECE", "-cpus=".$opt->{threads}, $opt->{proteinA}, $opt->{proteinB});
    run_cmd($L, \@cmd, $opt->{out});

    $opt->{proteinortho} = $opt->{basedir}."microPIECE.proteinortho";
}

sub run_cmd
{
    my ($log, $cmd, $path) = @_;

    my $currentdir = undef;

    if ($path)
    {
	$currentdir = getcwd;
	chdir($path);
    }

    $log->info("Calling command: ".join(" ", @{$cmd}));

    my $output = qx(@{$cmd});

    if ($? != 0)
    {
	$log->logdie("Error calling command: ".join(" ", @{$cmd}));
    }

    $log->debug("Output of command was: ".$output);

    if ($currentdir)
    {
	chdir($currentdir);
    }

    return $output;
}

=pod

Running the target prediction part of the microPIECE

=cut

sub run_targetprediction {

    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    unless (exists $opt->{run_targetprediction} && $opt->{run_targetprediction})
    {
	$L->info("Skipping target prediction step due to missing parameters");
	return;
    }

    $L->info("Starting target prediction step");

    $L->info("Finished target prediction step");

}


1;
