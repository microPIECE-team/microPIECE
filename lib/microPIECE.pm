package microPIECE;

use strict;
use warnings;

use version 0.77; our $VERSION = version->declare("v0.9.0");

use Log::Log4perl;
use Data::Dumper;
use Cwd;
use File::Path;
use File::Basename;
use File::Spec;
use File::Temp qw/ :POSIX /;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

use IPC::Run3;

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

sub check_files
{
    my $opt = shift;
    my @fileparameter = @_;

    foreach my $param (@fileparameter)
    {
	my $L = Log::Log4perl::get_logger();
	next unless (exists $opt->{$param} && defined $opt->{$param});
	if (-e $opt->{$param})
	{
	    $opt->{$param} = File::Spec->rel2abs($opt->{$param});
	} else {
	    $L->logdie(sprintf("Missing parameter for --%s or file '%s' is inaccessable\n", $param, $opt->{$param}));
	}

    }
}

=pod

Checking requirements for the pipeline

=cut

sub check_requirements {
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    $opt->{basedir} = getcwd."/";

    ##############################################################
    #
    #  Check settings for CLIP
    #
    ##############################################################

    # first we check if mandatory parameters have been specified
    check_files($opt, qw(genomeA genomeB annotationA annotationB));

    if (exists $opt->{clip} && @{$opt->{clip}} > 0)
    {
	for( my $i=0; $i<@{$opt->{clip}}; $i++ )
	{
	    my $file = $opt->{clip}[$i];
	    $L->logdie("Missing parameter for --clip or file '$file' is inaccessable\n") unless (-e $file);
	    $opt->{clip}[$i] = File::Spec->rel2abs($file);
	}
    } else {
	$L->logdie("Missing parameter for --clip\n");
    }

    unless (exists $opt->{adapterclip} && length($opt->{adapterclip})>0 && $opt->{adapterclip} =~ /^[ACGT]+$/i)
    {
	$L->logdie("Missing parameter for --adapterclip or unexpected characters provided");
    }

    # we need to run clip
    $opt->{run_clip} = 1;

    ##############################################################
    #
    #  Check settings for mining
    #
    ##############################################################

    # check if we should run mining
    if (exists $opt->{smallrnaseq} && (keys %{$opt->{smallrnaseq}}) > 0)
    {
	foreach my $condition (keys %{$opt->{smallrnaseq}})
	{
	    for( my $i=0; $i<@{$opt->{smallrnaseq}{$condition}}; $i++ )
	    {
		my $file = $opt->{smallrnaseq}{$condition}[$i];
		$L->logdie("Missing parameter for --smallrnaseq or file '$file' is inaccessable\n") unless (-e $file);
		$opt->{smallrnaseq}{$condition}[$i] = File::Spec->rel2abs($file);
	    }
	}
    }

    # we need to run mining
    $opt->{run_mining} = 1 if ((keys %{$opt->{smallrnaseq}}) > 0);

    # check for adapter sequences
    foreach my $adapter (qw(adaptersmallrnaseq5 adaptersmallrnaseq3))
    {
	unless (exists $opt->{$adapter} && defined ($opt->{$adapter}))
	{
	    if ($opt->{run_mining})
	    {
		$L->info(sprintf("No adapter sequence provided via --%s parameter, therefore this side will not be clipped\n", $adapter));
	    }
	    next;
	}
	if ((length($opt->{$adapter})>0) && $opt->{$adapter} =~ /[^ACGT]/i)
	{
	    $L->logdie(sprintf("Parameter --%s contains unexpected characters\n",$adapter));
	}
    }

    if ($opt->{run_mining})
    {
	# check for filter file
	check_files($opt, "filterncrnas");

	unless (exists  $opt->{speciesB_tag} &&
	        defined $opt->{speciesB_tag} &&
	        length( $opt->{speciesB_tag})>0)
	{
	    $L->logdie("Parameter --speciesB need to be specified for mining");
	}
    }

    ##############################################################
    #
    #  Check settings for target prediction
    #
    ##############################################################

    # check if we need to run the target prediction
    check_files($opt, "mirna");

    # if the mirna option is given, we will skip the mining process
    if ($opt->{mirna})
    {
	$L->info("Due to a provided miRNA database (--mirna parameter) we will skip the mining process");
	$opt->{run_mining} = 0;
    }

    # we need to run the targetprediction
    if (($opt->{run_clip} && $opt->{mirna}) || ($opt->{run_clip} && $opt->{run_mining}))
    {
	$opt->{run_targetprediction} = 1;
    }

    ##############################################################
    #
    #  Check output directory
    #
    ##############################################################

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

    ##############################################################
    #
    #  Switch into output directory
    #
    ##############################################################

    $opt->{basedir} = $opt->{basedir}.$opt->{out}."/";
    chdir($opt->{basedir});

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

    run_mining_clipping($opt);
    run_mining_filtering($opt);
    run_mining_downloads($opt);
    run_mining_mirbase_files($opt);
    run_mining_mirdeep2($opt);
    run_mining_complete($opt);
    run_mining_mirdeep2fasta($opt);
    run_mining_quantification($opt);

    run_mining_genomicposition($opt);

    $opt->{mirna} = $opt->{final_mature};
    $L->info("Finished mining step");
}

sub run_mining_genomicposition
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    # build a blast database
    my @cmd = ("makeblastdb", "-in", $opt->{genomeB}, "-dbtype", "nucl");
    run_cmd($L, \@cmd);

    # run the blast
    @cmd = ("blastn", "-query", $opt->{final_hairpin}, "-db", $opt->{genomeB}, "-num_threads", $opt->{threads}, "-dust", "no", "-soft_masking", "false", "-outfmt", '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore');
    my $blastoutput = run_cmd($L, \@cmd);

    # filter the hits and store them
    $opt->{mining}{genomic_location} = "miRNA_genomic_position.csv";
    open(FH, ">", $opt->{mining}{genomic_location}) || $L->logdie("Unable to open file '$opt->{mining}{genomic_location}': $!");
    foreach my $line (split(/\n/, $blastoutput))
    {
	my @fields = split(/\t/, $line);

	next unless ($fields[2] == 100                 # 100% identity
		     && $fields[3] == $fields[4]);     # hit over complete length
	print FH $line, "\n";
    }
    close(FH) || $L->logdie("Unable to close file '$opt->{mining}{genomic_location}': $!");


}

sub run_mining_quantification
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    # build an index for mapping
    my @cmd = ("bwa", "index", $opt->{final_mature});
    run_cmd($L, \@cmd);

    # map the filtered short read libraries
    my $temp_sai = tmpnam();
    my $temp_sam = tmpnam();
    my $temp_sam_mapped_only = tmpnam();
    my @config_file_content = ();

    foreach my $condition (keys %{$opt->{mining}{filtered}})
    {
	$opt->{mining}{quantification}{$condition} = [];
	foreach my $file (@{$opt->{mining}{filtered}{$condition}})
	{
	    my @cmd = ("bwa", "aln", "-n", 1, "-o", 0, "-e", 0, "-k", 1, "-t", $opt->{threads}, "-f", $temp_sai, $opt->{final_mature}, $file);
	    run_cmd($L, \@cmd);

	    # convert sai to sam
	    @cmd = ("bwa", "samse", "-f", $temp_sam, $opt->{final_mature}, $temp_sai, $file);
	    run_cmd($L, \@cmd);

	    # remove unmapped reads
	    @cmd = ("samtools", "view", "-F", 4, "-o", $temp_sam_mapped_only, $temp_sam);
	    run_cmd($L, \@cmd);

	    # expand multimappings
	    @cmd = ($opt->{scriptdir}."xa2multi.pl", $temp_sam_mapped_only);
	    my $final_sam_content = run_cmd($L, \@cmd);

	    # store the final output
	    my $final_filename = $file."_mapped.sam";
	    push(@{$opt->{mining}{quantification}{$condition}}, $final_filename);
	    push(@config_file_content, join("\t", $final_filename, $condition)."\n");

	    open(FH, ">", $final_filename) || $L->logdie("Unable to open file '$final_filename': $!");
	    print FH $final_sam_content;
	    close(FH) || $L->logdie("Unable to close file '$final_filename': $!");
	}
    }

    # write the config file
    my $config_file = "062_s2d_cfg";
    open(FH, ">", $config_file) || $L->logdie("Unable to open file '$config_file': $!");
    print FH @config_file_content;
    close(FH) || $L->logdie("Unable to close file '$config_file': $!");

    # run the quantification analysis
    $opt->{mining_quantification_result} = "miRNA_expression.csv";
    @cmd = ($opt->{scriptdir}."061_sam2de.pl", "--cfg", $config_file, "--mature_file", $opt->{final_mature}, "--out", $opt->{mining_quantification_result});
    run_cmd($L, \@cmd);
}

sub run_mining_mirdeep2fasta
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    $opt->{novel_mature} = "novel_mature.fa";
    $opt->{novel_hairpin} = "novel_hairpin.fa";

    my @cmd = ($opt->{scriptdir}."041_curated_mirdeep2fasta.pl", "--csv", $opt->{mirdeep_output}, "--cutoff", 10, "--matureout", $opt->{novel_mature}, "--hairpinout", $opt->{novel_hairpin}, "--species", $opt->{speciesB_tag});
    run_cmd($L, \@cmd);

    # combine novel and known mature sequences and ensure DNA nucleotides
    $opt->{final_mature} = "mature_combined_mirbase_novel.fa";
    open(OUT, ">", $opt->{final_mature}) || $L->logdie("Unable to open file '$opt->{final_mature}' for writing: $!");
    foreach my $file ($opt->{novel_mature}, "mature_mirbase.fa")
    {
	open(FH, "<", $file) || $L->logdie("Unable to open file '$file': $!");
	while(<FH>)
	{
	    if ($_ !~ /^>/)
	    {
		$_ =~ tr/Uu/Tt/;
	    }
	    print OUT $_;
	}
	close(FH) || $L->logdie("Unable to close file '$file': $!");
    }
    close(OUT) || $L->logdie("Unable to close file '$opt->{final_mature}' after writing: $!");

    # combine novel and known hairpin sequences and ensure DNA nucleotides
    $opt->{final_hairpin} = "hairpin_combined_mirbase_novel.fa";
    open(OUT, ">", $opt->{final_hairpin}) || $L->logdie("Unable to open file '$opt->{final_hairpin}' for writing: $!");
    foreach my $file ($opt->{novel_hairpin}, "precursor_mirbase.fa")
    {
	open(FH, "<", $file) || $L->logdie("Unable to open file '$file': $!");
	while(<FH>)
	{
	    if ($_ !~ /^>/)
	    {
		$_ =~ tr/Uu/Tt/;
	    }
	    print OUT $_;
	}
	close(FH) || $L->logdie("Unable to close file '$file': $!");
    }
    close(OUT) || $L->logdie("Unable to close file '$opt->{final_hairpin}' after writing: $!");
}

sub run_mining_complete
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my @cmd = ($opt->{scriptdir}."021_parse_miRDeep2_output.pl", "-mirdeep_out", $opt->{mirdeep_output}, "-mature_fasta", "mature_mirbase.fa");
    my $output = run_cmd($L, \@cmd);

    my $mirbase_completed = "mature_mirbase_completed.fa";
    open(FH, ">", $mirbase_completed) || $L->logdie("Unable to open file '$mirbase_completed': $!");
    print FH $output;
    close(FH) || $L->logdie("Unable to close file '$mirbase_completed': $!");
}

sub run_mining_rna2dna
{
    my ($inseq) = @_;

    my @lines = split(/\n/, ${$inseq});

    for(my $i=0; $i<@lines; $i++)
    {
	next if ($lines[$i] =~ /^>/);

	$lines[$i] =~ tr/Uu/Tt/;   # convert [Uu]->[Tt]
    }

    ${$inseq} = join("\n", @lines);
}

sub run_mining_mirdeep2
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();
    # miRDeep2 needs fasta headers without whitespaces and provides this script for that purpose

    my $genome_wo_whitespace = $opt->{basedir}."genome_without_whitespace.fa";
    my %files_from_mirbase = ();

    my @cmd=("remove_white_space_in_id.pl", $opt->{genomeB});
    my $genome_without_whitespace = run_cmd($L, \@cmd);
    # write it to its own file
    open(FH, ">", $genome_wo_whitespace) || $L->logdie("Unable to open file '$genome_wo_whitespace' for writing: $!");
    print FH $genome_without_whitespace;
    close(FH) || $L->logdie("Unable to close file '$genome_wo_whitespace' after writing: $!");
    # save some memory
    $genome_without_whitespace = "";

    # convert files from mirbase_files subroutine
    foreach my $file (qw(mature_mirbase.fa precursor_mirbase.fa mature.fa-no-speciesB.fa))
    {
	@cmd=("remove_white_space_in_id.pl", $file);
	my $output = run_cmd($L, \@cmd);
	# convert to DNA
	run_mining_rna2dna(\$output);
	# write it into a file
	my $new_filename = basename($file, ".fa")."_wo_whitespace.fa";
	open(FH, ">", $new_filename) || $L->logdie("Unable to open file '$new_filename' for writing: $!");
	print FH $output;
	close(FH) || $L->logdie("Unable to close file '$new_filename' after writing: $!");

	if (-s $new_filename > 0)
	{
	    # assume content
	    $files_from_mirbase{$file} = $new_filename;
	} else {
	    # assume not given and there none
	    $files_from_mirbase{$file} = "none";
	}
    }

    # create a bowtie index
    my $genome_index = "genome4mirdeep2";
    @cmd = ("bowtie-build", $genome_wo_whitespace, $genome_index);
    run_cmd($L, \@cmd);

    # run through all short read files and concat the output in a single fasta file
    my $tempfasta = tmpnam();
    foreach my $condition (keys %{$opt->{mining}{filtered}})
    {
	foreach my $file (@{$opt->{mining}{filtered}{$condition}})
	{
	    # convert to fasta
	    my @cmd = ("fastq2fasta.pl", $file);
	    my $output = run_cmd($L, \@cmd);
	    # write to tempfasta
	    open(FH, ">>", $tempfasta) || $L->logdie("Unable to open file '$tempfasta' for writing: $!");
	    print FH $output;
	    close(FH) || $L->logdie("Unable to close file '$tempfasta' after writing: $!");
	}
    }

    my $tempfasta_wo_whitespace = tmpnam();
    my $tempfasta_wo_whitespace_collapsed = tmpnam();
    my $temparf = tmpnam();

    # remove whitespaces
    @cmd=("remove_white_space_in_id.pl", $tempfasta);
    my $output = run_cmd($L, \@cmd);
    # write to tempfasta_wo_whitespace
    open(FH, ">", $tempfasta_wo_whitespace) || $L->logdie("Unable to open file '$tempfasta_wo_whitespace' for writing: $!");
    print FH $output;
    close(FH) || $L->logdie("Unable to close file '$tempfasta_wo_whitespace' after writing: $!");

    # collapse reads
    @cmd=("collapse_reads_md.pl", $tempfasta_wo_whitespace, $opt->{speciesB_tag});
    $output = run_cmd($L, \@cmd);
    # write to tempfasta_wo_whitespace_collapsed
    open(FH, ">", $tempfasta_wo_whitespace_collapsed) || $L->logdie("Unable to open file '$tempfasta_wo_whitespace_collapsed' for writing: $!");
    print FH $output;
    close(FH) || $L->logdie("Unable to close file '$tempfasta_wo_whitespace_collapsed' after writing: $!");

    # map to genome
    @cmd = ("mapper.pl", $tempfasta_wo_whitespace_collapsed, "-c", "-q", "-n", "-l", 17, "-p", $genome_index, "-t", $temparf);
    run_cmd($L, \@cmd);

    # run mirdeep
    @cmd = ("miRDeep2.pl", $tempfasta_wo_whitespace_collapsed, $genome_wo_whitespace, $temparf, $files_from_mirbase{'mature_mirbase.fa'}, $files_from_mirbase{'mature.fa-no-speciesB.fa'}, $files_from_mirbase{'precursor_mirbase.fa'}, "-P");
    run_cmd($L, \@cmd);

    # save output
    my @csv_files = grep { /result_/ } (glob("*.csv"));
    unless (@csv_files == 1)
    {
	$L->logdie("Found multiple *.csv files instead of a single: ".join(", ", map {"'$_'"} (@csv_files)));
    }

    $opt->{mirdeep_output} = "mirdeep_output.csv";
    rename $csv_files[0], $opt->{mirdeep_output} || $L->logdie("Unable to rename mirdeep output file");
}

sub run_mining_mirbase_files
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();
    # -species	 	:= 3letter code of desired species(b)
    # -precursor_file 	:= hairpin.fasta from miRBase.org
    # -mature 		:= mature.fasta from miRBase.org
    # -organism 	:= organism.txt from miRBase.org
    # -out 		:= output folder
    # This script separates the miRBase files into groups of multifasta files that either belong to the speciesB or not.
    # In every case, it filters the microRNAs so that only metazoan are included.
    my @cmd = ($opt->{scriptdir}."011_mirbase_files.pl",
	       "-species", $opt->{speciesB_tag},
	       "-precursor_file", $opt->{mining}{download}{hairpin},
	       "-mature", $opt->{mining}{download}{mature},
	       "-organism", $opt->{mining}{download}{organisms},
	       "-out", "./");
    run_cmd($L, \@cmd);
}

sub run_mining_downloads
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my %filelist = (
	organisms => "organisms.txt.gz",
	mature    => "mature.fa.gz",
	hairpin   => "hairpin.fa.gz"
	);

    foreach my $key(sort keys %filelist)
    {
	my $file = $filelist{$key};
	my @cmd=("wget", "--quiet", "ftp://mirbase.org/pub/mirbase/CURRENT/".$file);
	run_cmd($L, \@cmd);
	# decompress the file
	my $output = basename($file, ".gz");
	my $status = gunzip $file => $output || $L->logdie("gunzip failed: $GunzipError");

	$opt->{mining}{download}{$key} = $output;
    }
}

# If a set of ncRNAs (with exception of miRNAs, of course) is provided,
# the smallRNA reads are filtered against those ncRNAs to exclude unwanted fragments.
sub run_mining_filtering
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    unless ($opt->{filterncrnas})
    {
	# nothing to do, therefore skip the step
	$L->info("Filtering step was skipped, due to missing --filterncrnas parameter");
	$opt->{mining}{filtered} = $opt->{mining}{trimmed};
	return;
    }

    my @cmd = ("bwa", "index", $opt->{filterncrnas});
    run_cmd($L, \@cmd);

    $L->logdie("Something went completly wrong") unless (exists $opt->{mining} && $opt->{mining}{trimmed});
    my $tempbwaout = tmpnam();
    my $tempsamout = tmpnam();
    my $tempsamfilteredout = tmpnam();
    my $tempsamsortedout = tmpnam();
    foreach my $condition (keys %{$opt->{mining}{trimmed}})
    {
	foreach my $file (@{$opt->{mining}{trimmed}{$condition}})
	{
	    # new filename will be
	    my $filtered_fq = $opt->{basedir}.basename($file, (".fq", ".fastq"))."_filtered.fq";

	    # map to filter file
	    # -n 1 := edit distance of 1
	    # -o 0 := no gap opens
	    # -e 0 := no gap extensions
	    # -k 1 := max 1 difference in the seed
	    my @cmd = ("bwa", "aln", "-n", 1, "-o", 0,
		       "-e", 0, "-k", 1, "-t", $opt->{threads},
		       "-f", $tempbwaout, $opt->{filterncrnas}, $file);
	    run_cmd($L, \@cmd);
	    # convert to sam
	    @cmd = ("bwa", "samse", "-f", $tempsamout, $opt->{filterncrnas}, $tempbwaout, $file);
	    run_cmd($L, \@cmd);
	    # filter sam
	    @cmd = ("samtools", "view",
		    "--threads", $opt->{threads},
		    "-b", "-o", $tempsamfilteredout, "-f", 4, $tempsamout);
	    run_cmd($L, \@cmd);
	    # sort sam
	    @cmd = ("samtools", "sort",
		    "--threads", $opt->{threads},
		    "-o", $tempsamsortedout, $tempsamfilteredout);
	    run_cmd($L, \@cmd);
	    # extract fastqs
	    @cmd = ("bedtools", "bamtofastq", "-i", $tempsamsortedout, "-fq", $filtered_fq);
	    run_cmd($L, \@cmd);

	    push(@{$opt->{mining}{filtered}{$condition}}, $filtered_fq);
	}
    }

    return;
}

sub run_mining_clipping
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my $minlen = 17;
    # --trim-n		:= trim N's at the read ends
    # --minimum-length	:= minimal read length of 17 nt
    my @cmd = ("cutadapt", "--trim-n", "--minimum-length", $minlen);
    if ($opt->{adaptersmallrnaseq3})
    {
	push(@cmd, ("--adapter", $opt->{adaptersmallrnaseq3}))
    }
    if ($opt->{adaptersmallrnaseq5})
    {
	push(@cmd, ("--front", $opt->{adaptersmallrnaseq5}))
    }

    foreach my $condition (keys %{$opt->{smallrnaseq}})
    {
	foreach my $file (@{$opt->{smallrnaseq}{$condition}})
	{
	    my @cmd2run = @cmd;
	    my $outfile = $opt->{basedir}.basename($file, (".fq", ".fastq"))."_trimmed.fq";
	    push(@cmd2run, ($file, "-o", $outfile));
	    run_cmd($L, \@cmd2run);
	    push(@{$opt->{mining}{trimmed}{$condition}}, $outfile);
	}
    }
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

    foreach my $file (@{$opt->{seq4prediction}})
    {
	my $final_output = $opt->{basedir}.basename($file)."_final_miranda_output.txt";

	my @cmd = ($opt->{scriptdir}."Targetprediction.pl", $opt->{mirna}, $file, $final_output);
	run_cmd($L, \@cmd)
    }
    $L->info("Finished target prediction step");

}

sub run_CLIP_transfer
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    $opt->{mRNAB} = $opt->{basedir}."/"."mRNAB.fa";
    my @cmd = ("gffread", $opt->{annotationB}, "-w", $opt->{mRNAB}, "-F", "-g", $opt->{genomeB});
    run_cmd($L, \@cmd);

    my @inputfiles = glob("clip_merged_*of*BEDfilter_mapGFF_minLen*_min*_max*_sort_UC.fasta");

    my $file_uniqueA     = $opt->{basedir}.basename($opt->{annotationA}, ".gff")."_unique.csv";
    my $file_uniqueA_log = $opt->{basedir}.basename($opt->{annotationA}, ".gff")."_unique.err";
    @cmd = ($opt->{scriptdir}."CLIP_parse_gff_return_longest_transcript.pl", $opt->{annotationA});
    run_cmd($L, \@cmd, undef, $file_uniqueA, $file_uniqueA_log);

    my $file_uniqueB     = $opt->{basedir}.basename($opt->{annotationB}, ".gff")."_unique.csv";
    my $file_uniqueB_log = $opt->{basedir}.basename($opt->{annotationB}, ".gff")."_unique.err";
    @cmd = ($opt->{scriptdir}."CLIP_parse_gff_return_longest_transcript.pl", $opt->{annotationB});
    run_cmd($L, \@cmd, undef, $file_uniqueB, $file_uniqueB_log);

    foreach my $file (@inputfiles)
    {
	my $needle_csv = $opt->{basedir}.basename($file, ".fasta")."_needle.csv";
	my $needle_aln = $opt->{basedir}.basename($file, ".fasta")."_needle.aln";
	my $bed_out = $opt->{basedir}.basename($file, ".fasta")."_transfered.bed";
	my $bed_merged = $opt->{basedir}.basename($file, ".fasta")."_transfered_merged.bed";
	my $final_fasta = $opt->{basedir}.basename($file, ".fasta")."_transfered_final.fasta";

	@cmd = ($opt->{scriptdir}."CLIP_map_clip_gff_needle.pl", $file_uniqueA, $opt->{proteinortho}, $file, $file_uniqueB, $opt->{mRNAB}, $needle_csv);
	run_cmd($L, \@cmd, undef, $needle_aln);

	# convert csv into bed file
	@cmd = ($opt->{scriptdir}."CLIP_csv_to_bed.pl", $needle_csv, $bed_out);
	run_cmd($L, \@cmd);

	# merge bed annotations
	@cmd = ("bedtools", "merge", "-i", $bed_out, "-c", 4, "-o", "collapse");
	run_cmd($L, \@cmd, undef, $bed_merged);

	# convert merged bed into sequences based on transcripts
	@cmd = ("bedtools", "getfasta", "-name", "-fi", $opt->{mRNAB}, "-bed", $bed_merged, "-fo", $final_fasta);
	run_cmd($L, \@cmd);

	push(@{$opt->{seq4prediction}}, $final_fasta);
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
	my @cmd = ($opt->{scriptdir}."CLIP_bedtool_discard_sizes.pl", $file, $min, $max);
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
	@cmd = ($opt->{scriptdir}."CLIP_fasta_uc_and_filter4annotations.pl", $fasta);
	run_cmd($L, \@cmd, undef, $fastaUC);
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
	my @cmd = ($opt->{scriptdir}."CLIP_mapper.pl", $file, $opt->{annotationA}, $minlength);
	run_cmd($L, \@cmd, undef, $outputname);
    }
}

sub run_CLIP_filterbed
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my @cmd = ($opt->{scriptdir}."CLIP_bed2signal.pl", $opt->{basedir}."clip_merged.bed");

    my $num_fields = int(@{$opt->{clip}});

    for(my $i=1;$i<=@{$opt->{clip}};$i++)
    {
	my $filtered_bed_out = sprintf("%sclip_merged_%dof%dBEDfilter.bed", $opt->{basedir}, $i, $num_fields);
	run_cmd($L, [@cmd, $i], undef, $filtered_bed_out);
    }
}


sub run_CLIP_bedtools_merge
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();

    my @cmd = ($opt->{scriptdir}."CLIP_merge_bed_files.pl", "--output", $opt->{basedir}."clip_merged.bed", "--log", "merging_bed_files.log");

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

    my $gsnap_output = tmpnam();
    my $samtools_output = tmpnam();

    foreach my $clipfile (@{$opt->{clip}})
    {
	my $trimmedfile = $opt->{basedir}.basename($clipfile).".trim";
	my $bamfile     = $opt->{basedir}.basename($clipfile).".bam";
	my $bedfile     = $opt->{basedir}.basename($clipfile).".bed";
	# -N 1		:= look for splice sites
	# -B 5		:= batch mode 5, allocate positions, genome and suffix array
	# -O		:= ordered output
	# -A sam 	:= output format

	my @cmd = ("gsnap", "-N", 1, "-B", 5, "-O", "-A", "sam", "-t", $opt->{threads}, "-D", "speciesA_db", "-d", "speciesA", $trimmedfile);
	run_cmd($L, \@cmd, undef, $gsnap_output);

	@cmd = ("samtools", "view", "-Sb", $gsnap_output);
	run_cmd($L, \@cmd, undef, $samtools_output);

	@cmd = ("samtools", "sort", "-o", $bamfile, $samtools_output);
	run_cmd($L, \@cmd);

	@cmd = ("bedtools", "bamtobed", "-i", $bamfile);
	run_cmd($L, \@cmd, undef, $bedfile);
    }
}

sub run_CLIP_build_db
{
    my ($opt) = @_;

    my $L = Log::Log4perl::get_logger();
    # -k 15 := kmer size 15 for genome index
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
	# -m 20		:= min length of read
	# --trim-n	:= trim terminal Ns of reads
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
    my ($log, $cmd, $path, $outfile, $errfile) = @_;

    my $currentdir = undef;

    if ($path)
    {
	$currentdir = getcwd;
	chdir($path);
    }

    my $redirection = "";

    my $stdout = "";
    my $stderr = "";
    my $param_stdout = \$stdout;
    if (defined $outfile)
    {
	$param_stdout = $outfile;
	$redirection .= " >$outfile";
    }
    my $param_stderr = \$stderr;
    if (defined $errfile)
    {
	$param_stderr = $errfile;
	$redirection .= " 2>$errfile";
    }

    $log->info("Calling command: ".join(" ", @{$cmd}, $redirection));

    run3($cmd, undef, $param_stdout, $param_stderr);
    #$stdout = qx(@{$cmd});

    if ($? != 0)
    {
	$log->logdie("Error calling command: ".join(" ", @{$cmd})."STDERR: $stderr");
    }

    $log->debug("Output of command was: ".$stdout);
    $log->debug("Output of command was: ".$stderr);

    if ($currentdir)
    {
	chdir($currentdir);
    }

    return $stdout;
}


1;
