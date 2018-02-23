package microPIECE;

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

    $L->info("Finished CLIP step");

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
