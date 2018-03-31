#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin $Script);
use lib "$RealBin/lib/";

use IO::Handle;

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
    adaptersmallrnaseq5=> undef,
    adaptersmallrnaseq3=> undef,
    filterncrnas       => undef,
    config             => undef,
    threads            => 1,
    out                => "out",
    overwrite          => 0,
    scriptdir          => $RealBin."/scripts/",
    testrun            => undef,
    mirna              => undef,
    speciesB_tag       => undef,
    mirbasedir         => undef,
    tempdir            => undef,
    piranha_bin_size   => 20,
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
    'adaptersmallrnaseq5=s'=> \$opt->{adaptersmallrnaseq5},
    'adaptersmallrnaseq3=s'=> \$opt->{adaptersmallrnaseq3},
    'filterncrnas=s'       => \$opt->{filterncrnas},
    'config=s'             => \$opt->{config},
    'threads=i'            => \$opt->{threads},
    'overwrite'            => \$opt->{overwrite},
    'testrun'              => \$opt->{testrun},
    'out=s'                => \$opt->{out},
    'mirna=s'              => \$opt->{mirna},
    'speciesBtag=s'        => \$opt->{speciesB_tag},
    'mirbasedir=s'         => \$opt->{mirbasedir},
    'tempdir=s'            => \$opt->{tempdir},
    'piranhabinsize=i'     => \$opt->{piranha_bin_size},
    ) || pod2usage(1);

# split clip files if required
$opt->{clip} = [ split(",", join(",", @{$opt->{clip}})) ];
# split rnaseq files if required
foreach my $cond (keys %{$opt->{smallrnaseq}})
{
    $opt->{smallrnaseq}{$cond} = [ split(",", join(",", @{$opt->{smallrnaseq}{$cond}})) ];
}

# if no tempdir was set, we will use $opt->{out}/tmp
unless (exists $opt->{tempdir} && defined $opt->{tempdir})
{
    $opt->{tempdir} = $opt->{out}."/tmp/";
}

# help
$opt->{help} && pod2usage(1);

# version
if($opt->{version}){
	print $microPIECE::VERSION->normal(),"\n";
	exit 0;
}

microPIECE::hello();
STDOUT->flush();

microPIECE::check_requirements($opt);

microPIECE::print_settings($opt);

microPIECE::run_mining($opt);

microPIECE::run_clip($opt);

microPIECE::run_targetprediction($opt);

microPIECE::transfer_resultfiles($opt);

__END__

=pod

=encoding utf8

=head1 NAME

microPIECE - microRNA pipeline enhanced by CLIP experiments

=head1 SYNOPSIS

 microPIECE.pl \
   --genomeA <speciesA_genome.fa> \
   --genomeB <speciesB_genome.fa> \
   --annotationA <speciesA_genome.gff> \
   --annotationB <speciesB_genome.gff> \
   --clip <ago_clip_seq.fastq> \
   --mirna <mature_miRNAs.fa>

=head1 DESCRIPTION

The microPIECE (microRNA pipeline enhanced by CLIP experiments)
takes the AGO-CLIP data from a speciesA and transfers it to a speciesB.
Given a set of miRNAs from speciesB it then predicts their targets on the transfered CLIP regions.

For the minimal workflow it needs a genome file, as well as its annotation file in GFF format for speciesA and speciesB.
For speciesA at least one AGO-CLIP dataset is needed and speciesB needs a set of miRNAs for the target prediction.
For the full workflow, a set of smallRNA-sequencing data is additionally needed and a set of non-coding RNAs can be
provided as filter. The pipeline uses the smallRNA data for the mining of novel microRNAs and the completion of
the given miRNA dataset, if needed. It further performs expression calculation, isoform detection, genomic loci
identification and orthology determination.

=head1 EXAMPLE

 microPIECE.pl \
     --genomeA testset/NC_035109.1_reduced_AAE_genome.fa \
     --genomeB testset/NC_007416.3_reduced_TCA_genome.fa \
     --annotationA testset/NC_035109.1_reduced_AAE_genome.gff \
     --annotationB testset/NC_007416.3_reduced_TCA_genome.gff \
     --clip testset/SRR5163632_aae_clip_reduced.fastq,testset/SRR5163633_aae_clip_reduced.fastq,testset/SRR5163634_aae_clip_reduced.fastq \
     --clip testset/SRR5163635_aae_clip_reduced.fastq,testset/SRR5163636_aae_clip_reduced.fastq,testset/SRR5163637_aae_clip_reduced.fastq \
     --adapterclip GTGTCAGTCACTTCCAGCGG \
     --overwrite \
     --smallrnaseq a=testset/tca_smallRNAseq_rna_contaminated.fastq \
     --adaptersmallrnaseq3=TGGAATTCTCGGGTGCCAAGG \
     --adaptersmallrnaseq5 GTTCAGAGTTCTACAGTCCGACGATC \
     --filterncrnas testset/TCA_all_ncRNA_but_miR.fa \
     --speciesB tca 2>&1 | tee out.log

=head1 PARAMETERS

=over 4

=item C<--version|-V>

version of this pipeline

=item C<--help|-h>

prints a helpful help message

=item C<--genomeA> and C<--genomeB>

Genome of the species with the CLIP data (species A, C<--genomeA>) and
the genome of the species where we want to predict the miRNA targets
(species B, C<--genomeB>)

=item C<--gffA> and C<--gffB>

Genome feature file (GFF) of the species with the CLIP data (species
A, C<--gffA>) and the GFF of the species where we want to predict the
miRNA targets (species B, C<--gffB>)

=item C<--clip>

Comma-separated CLIP-seq .fastq files in Format

  --clip con1_rep1_clip.fq,con1_rep2_clip.fq,con2_clip.fq
  # OR
  --clip con1_rep1_clip.fq --clip con1_rep2_clip.fq --clip con2_clip.fq

=item 	C<--adapterclip>

Sequencing-adapter of CLIP reads

=item C<--smallrnaseq>

Comma-separated smallRNA-seq FASTQ files, initialized with
'condition=' in Format

  --smallrnaseq con1=A.fastq,B.fastq --smallrnaseq con2=C.fq
  # OR
  --smallrnaseq con1=A.fastq --smallrnaseq con1=B.fastq --smallrnaseq con2=C.fq

=item C<--adaptersmallrnaseq5> and C<--adaptersmallrnaseq3>

5' adapter of smallRNA-seq reads (C<--adaptersmallrnaseq5>) and for 3' end (C<--adaptersmallrnaseq3>)

=item C<--filterncrnas>

Multi-fasta file of ncRNAs to filter smallRNA-seq reads. Those must
not contain miRNAs.

=item C<--threads>

Number of threads to be used

=item C<--overwrite>

set this parameter to overwrite existing files

=item C<--testrun>

sets this pipeline to testmode (accounting for small testset in
piranha). This option should not be used in real analysis!

=item C<--out>

output folder

=item C<--mirna>

miRNA set, if set, mining is disabled and this set is used for prediction

=item C<--speciesBtag>

Three letter code of species where we want to predict the miRNA
targets (species B, C<--speciesBtag>).

=item C<--mirbasedir>

The folder specified by C<--mirbasedir> is searched for the files
F<organisms.txt.gz>, F<mature.fa.gz>, and F<hairpin.fa.gz>. If the
files are not exist, they will be downloaded.

=item C<--tempdir>

The folder specified by C<--tempdir> is used for temporary files. The
default value is F<tmp/> inside the output folder specified by the
C<--out> parameter.

=item C<--piranhabinsize>

Sets the F<Piranha> bin size and has a default value of C<20>.

=back

=head1 OUTPUT

=over 4

=item pseudo mirBASE dat file: F<final_mirbase_pseudofile.dat>

A pseudo mirBASE dat file containing all precursor sequences with their named mature sequences and their coordinates. It only contain the fields:

=over 4

=item C<ID>

=item C<FH> and C<FT>

=item C<SQ>

=back

=item mature miRNA set: F<mature_combined_mirbase_novel.fa>

mature microRNA set, containing novels and miRBase-completed (if mined), together with the known miRNAs from miRBase

=item precursor miRNA set: F<hairpin_combined_mirbase_novel.fa>

precursor microRNA set, containing novels (if mined), together with the known miRNAs from miRBase

=item mature miRNA expression per condition: F<miRNA_expression.csv>

Semicolon-separated file containing:

=over 8

=item 1. C<rpm>

=item 2. C<condition>

=item 3. C<miRNA>

=back

=item orthologous prediction file: F<miRNA_orthologs.csv>

tab-separated file containing:

=over 8

=item 1. C<query_id>

=item 2. C<subject_id>

=item 3. C<identity>

=item 4. C<alignment length>

=item 5. C<number mismatches>

=item 6. C<number gap openings>

=item 7. C<start position inside query>

=item 8. C<end position inside query>

=item 9. C<start position inside subject>

=item 10. C<end position inside subject>

=item 11. C<evalue>

=item 12. C<bitscore>

=item 13. C<aligned query sequence>

=item 14. C<aligned subject sequence>

=item 15. C<length query sequence>

=item 16. C<length subject sequence>

=item 17. C<coverage for query sequence>

=item 18. C<coverage for subject sequence>

=back

=item miRDeep2 mining result in HTML/CSV F<mirdeep_output.html/csv>

the standard output HTML/CSV file of miRDeep2

=item ISOMIR prediction files: F<isomir_output_CONDITION.csv>

semincolon delimited file containing:

=over 8

=item 1. C<mirna>

=item 2. C<substitutions>

=item 3. C<added nucleotids on 3' end>

=item 4. C<nucleotides at 5' end different from the annonated sequence>

=item 5. C<nucleotides at 3' end different from the annonated sequence>

=item 6. C<sequence>

=item 7. C<rpm>

=item 8. C<condition>

=back

=item genomics location of miRNAs: F<miRNA_genomic_position.csv>

tab delimited file containing:

=over 8

=item 1. C<miRNA>

=item 2. C<genomic contig>

=item 3. C<identify>

=item 4. C<length>

=item 5. C<miRNA-length>

=item 6. C<number mismatches>

=item 7. C<number gapopens>

=item 8. C<miRNA-start>

=item 9. C<miRNA-stop>

=item 10. C<genomic-start>

=item 11. C<genomic-stop>

=item 12. C<evalue>

=item 13. C<bitscore>

=back

=item all library support-level target predictions: F<*_miranda_output.txt>

miranda output, reduced to the lines, starting with > only

=item all library support-level CLIP transfer .bed files: F<*transfered_merged.bed>

bed-file of the transferred CLIP-regions in speciesB transcriptome

=back

=head1 CAVEATS

Complete list of open issues is available on L<Github-Issues|https://github.com/microPIECE-team/microPIECE/issues>.

Please report any new issues ad L<new Github-Issue|https://github.com/microPIECE-team/microPIECE/issues/new>.

=head1 CHANGELOG

=over 4

=item scheduled for next release

No features planed

=item L<v1.4.0|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.4.0> (2018-03-31)

Copying pseudo mirBASE dat file C<final_mirbase_pseudofile.dat> into output folder (Fixes L<#131|https://github.com/microPIECE-team/microPIECE/issues/131>)

Corrected C<RNA::HairpinFigure> output (Fixes L<#137|https://github.com/microPIECE-team/microPIECE/issues/137>)

Fix the requirement of an accession inside mirBASE dat file (Fixes L<#134|https://github.com/microPIECE-team/microPIECE/issues/134>)

Avoiding error message while copying the out file for genomic location into base folder (Fixes L<#117|https://github.com/microPIECE-team/microPIECE/issues/117>)

=item L<v1.3.0|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.3.0> (2018-03-29)

Creating all structures on the fly using pseudo-mirBASE-dat as input.

Using F<miRNA.dat> from mirBASE as source for mature/precursor sequence and relationship (Fixes L<#127|https://github.com/microPIECE-team/microPIECE/issues/127>)

Fix of division-by-zero bug for empty mapping files (Fixes L<#118|https://github.com/microPIECE-team/microPIECE/issues/118>)

Fix of typo in C<--piranhabinsize> option (Fixes L<#116|https://github.com/microPIECE-team/microPIECE/issues/116>)

=item L<v1.2.3|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.3> (2018-03-26)

Fix transformation of precursor sequences based on mirbase #22 precursor sequences with a single mature.
(Fixes L<#109|https://github.com/microPIECE-team/microPIECE/issues/109>)

=item L<v1.2.2|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.2> (2018-03-23)

Improved collision detection for newly identified miRNAs avoiding crashed caused by genomic copies.
(Fixes L<#105|https://github.com/microPIECE-team/microPIECE/issues/105>)

=item L<v1.2.1|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.1> (2018-03-23)

Enables stable numbering for newly identified miRNAs based on their precursor and mature sequences
(Fixes L<#101|https://github.com/microPIECE-team/microPIECE/issues/101>)

=item L<v1.2.0|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.0> (2018-03-22)

We are using miraligner which requires a java version 1.7, but 1.8 was
installed by default. This was fixed by switching to v1.4 of the
docker base image. Additionally, miraligner requires fix filenames for
its databases. Therefore, the version v1.2.0 solved miraligner related
bugs and reenables the isomir detection.  (Fixes
L<#97|https://github.com/microPIECE-team/microPIECE/issues/97> and
L<#98|https://github.com/microPIECE-team/microPIECE/issues/98>)

=item L<v1.1.0|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.1.0> (2018-03-12)

Add isomir detection and copy the final genomic location file to the
output filter (Fixes
L<#34|https://github.com/microPIECE-team/microPIECE/issues/34>)

=item L<v1.0.7|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.7> (2018-03-08)

Piranha was lacking of a bin_size parameter. Added parameter C<--piranhabinsize> with a default value of C<20>
(Fixes L<#66|https://github.com/microPIECE-team/microPIECE/issues/80>)

=item L<v1.0.6|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.6> (2018-03-08)

Added parameter C<--mirbasedir> and C<--tempdir> to support local
mirbase files and relocation of directory for temporary files (Fixes
L<#66|https://github.com/microPIECE-team/microPIECE/issues/66>,
L<#73|https://github.com/microPIECE-team/microPIECE/issues/73>, and
L<#76|https://github.com/microPIECE-team/microPIECE/issues/76>)

=item L<v1.0.5|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.5> (2018-03-07)

Update of documentation and correct spelling of C<--mirna> parameter

=item L<v1.0.4|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.4> (2018-03-07)

Fixes complete mature in final output (Fixes L<#69|https://github.com/microPIECE-team/microPIECE/issues/69>)

=item L<v1.0.3|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.3> (2018-03-06)

Add tests for perl scripts in script folder which ensure the correct handling of BED stop coordinates (Fixes L<#65|https://github.com/microPIECE-team/microPIECE/issues/65>)

=item L<v1.0.2|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.2> (2018-03-05)

Fixes the incorrect sorting of BED files, result was correct, but sorting was performed in the wrong order. (Fixes L<#63|https://github.com/microPIECE-team/microPIECE/issues/63>)

=item L<v1.0.1|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.1> (2018-03-05)

Fix an error conserning BED file handling of start and stop coordinates. (Fixes L<#59|https://github.com/microPIECE-team/microPIECE/issues/59>)

=item L<v1.0.0|https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.0> (2018-03-05)

=begin html

is archived as <a href="https://doi.org/10.5281/zenodo.1188484"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1188484.svg" alt="DOI" /></a>
and submitted to <a href="http://joss.theoj.org">The Journal of Open Source Software</a>.

=end html

=begin text

is archived as L<https://doi.org/10.5281/zenodo.1188484> and submitted to L<The Journal of Open Source Software|https://joss.theoj.org>.

=end text

=item L<v0.9.0|https://github.com/microPIECE-team/microPIECE/releases/tag/v0.9.0> (2018-03-05)

=begin html

first version archived at Zenodo with the <a href="https://doi.org/10.5281/zenodo.1188481"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1188481.svg" alt="DOI" /></a>

=end html

=begin text

first version archived at Zenodo with the  L<https://doi.org/10.5281/zenodo.1188481>

=end text

=back

=head1 LICENCE

This program is released under GPLv2.
For further license information, see F<LICENSE.md> shipped with this program.
Copyright(c)2018 Daniel Amsel and Frank Förster (employees of Fraunhofer Institute for Molecular Biology and Applied Ecology IME).
All rights reserved.

=head1 AUTHORS

=over 4

=item Daniel Amsel E<lt>daniel.amsel@ime.fraunhofer.deE<gt>

=item Frank Förster E<lt>frank.foerster@ime.fraunhofer.deE<gt>

=back

=head1 SEE ALSO

L<Project source code on Github|https://github.com/microPIECE-team/microPIECE>
L<Docker image on DockerHub|https://hub.docker.com/r/micropiece/micropiece/>
L<Travis continuous integration page|https://travis-ci.org/microPIECE-team/microPIECE>
L<Test coverage reports|https://coveralls.io/github/microPIECE-team/microPIECE>

=cut
