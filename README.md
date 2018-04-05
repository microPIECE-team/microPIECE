# microPIECE

The `microPIECE` (**micro**RNA **pi**peline **e**nhanced by **C**LIP **e**xperiments) takes the AGO-CLIP data from a *speciesA* and transfers it to a *speciesB*. Given a set of miRNAs from *speciesB* it then predicts their targets on the transfered CLIP regions.

For the minimal workflow it needs a genome file, as well as its annotation file in GFF format for *speciesA* and *speciesB*. For *speciesA* at least one AGO-CLIP dataset is needed and *speciesB* needs a set of miRNAs for the target prediction. For the full workflow, a set of smallRNA-sequencing data is additionally needed and a set of non-coding RNAs can be provided as filter. The pipeline uses the smallRNA data for the mining of novel microRNAs and the completion of the given miRNA dataset, if needed. It further performs expression calculation, isoform detection, genomic loci identification and orthology determination.
## Status
[![Build Status](https://travis-ci.org/microPIECE-team/microPIECE.svg?branch=master)](https://travis-ci.org/microPIECE-team/microPIECE)
[![Coverage Status](https://coveralls.io/repos/github/microPIECE-team/microPIECE/badge.svg?branch=travis)](https://coveralls.io/github/microPIECE-team/microPIECE?branch=travis)

## Required Software
  - [bwa](http://bio-bwa.sourceforge.net/) (0.7.12-r1039)
  - [samtools](http://samtools.sourceforge.net/) (1.4.1)
  - [bedtools](http://bedtools.readthedocs.io/en/latest/) (2.27.1)
  - [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (1.1.2)
  - [miRDeep2](https://www.mdc-berlin.de/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/documentation) (2.0.0.8)
  - [miraligner](https://github.com/lpantano/seqcluster) (1.2.4a)
  - [NCBI-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (2.2.31+)
  - [Proteinortho](https://www.bioinf.uni-leipzig.de/Software/proteinortho/) (5.16b) 
  - [Cutadapt](https://github.com/marcelm/cutadapt) (1.9.1)
  - [gmap/gsnap](http://research-pub.gene.com/gmap/) (2018-02-12)
  - [Piranha](http://smithlabresearch.org/software/piranha/) (1.2.1)
  - [miranda](http://34.236.212.39/microrna/getDownloads.do) (aug2010)

## Required Perl modules
  - [Getopt::Long](http://search.cpan.org/dist/Getopt-Long/lib/Getopt/Long.pm) (2.5)
  - [File::Temp](http://search.cpan.org/~dagolden/File-Temp-0.2304/lib/File/Temp.pm) (0.2304)
  - [RNA::HairpinFigure](http://search.cpan.org/~shenwei/RNA-HairpinFigure-0.141212/lib/RNA/HairpinFigure.pm) (0.141212)
  - [Pod::Usage](http://search.cpan.org/~marekr/Pod-Usage-1.69/lib/Pod/Usage.pm) (1.69)
  - [Log::Log4perl](http://search.cpan.org/~mschilli/Log-Log4perl-1.49/lib/Log/Log4perl.pm) (1.49)

## Installation
Please install the dependencies and run
```
git clone -b v1.4.0 git@github.com:microPIECE-team/microPIECE.git
```
or download the latest release as `*.tar.gz` or `*.zip` file:
```
curl -L -o microPIECE_v1.4.0.tar.gz https://github.com/microPIECE-team/microPIECE/archive/v1.4.0.tar.gz
# or
curl -L -o microPIECE_v1.4.0.zip https://github.com/microPIECE-team/microPIECE/archive/v1.4.0.zip
```

## Docker
We also provide `microPIECE` as [DOCKER image](https://hub.docker.com/r/micropiece/micropiece/). We tested the image on Ubuntu, Debian and MacOS. For the latter one, the `Piranha` command `make test` fails during the build, but when entering the container, the test succeds. Therefore, we temporarily excluded this statement.

### Information about the docker images:
| Branch | Size | Layers | Comment |
|-|-|-|-|
|[![](https://images.microbadger.com/badges/version/micropiece/micropiece:v1.4.0.svg)](https://microbadger.com/images/micropiece/micropiece:v1.4.0) | [![](https://images.microbadger.com/badges/image/micropiece/micropiece:v1.4.0.svg)](https://microbadger.com/images/micropiece/micropiece:v1.4.0) | [![](https://images.microbadger.com/badges/commit/micropiece/micropiece:v1.4.0.svg)](https://microbadger.com/images/micropiece/micropiece:v1.4.0) | Latest release [v1.4.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.4.0) |
|[![](https://images.microbadger.com/badges/version/micropiece/micropiece:master.svg)](https://microbadger.com/images/micropiece/micropiece:master) | [![](https://images.microbadger.com/badges/image/micropiece/micropiece:master.svg)](https://microbadger.com/images/micropiece/micropiece:master) | [![](https://images.microbadger.com/badges/commit/micropiece/micropiece:master.svg)](https://microbadger.com/images/micropiece/micropiece:master) | |
|[![](https://images.microbadger.com/badges/version/micropiece/micropiece:develop.svg)](https://microbadger.com/images/micropiece/micropiece:develop) | [![](https://images.microbadger.com/badges/image/micropiece/micropiece:develop.svg)](https://microbadger.com/images/micropiece/micropiece:develop) | [![](https://images.microbadger.com/badges/commit/micropiece/micropiece:develop.svg)](https://microbadger.com/images/micropiece/micropiece:develop) | |


```
docker pull micropiece/micropiece:v1.4.0
git clone git@github.com:microPIECE-team/microPIECE-testset.git testset
docker run -it --rm -v $PWD:/data micropiece/micropiece:v1.4.0 microPIECE.pl   \
  --genomeA testset/NC_035109.1_reduced_AAE_genome.fa  \
  --genomeB testset/NC_007416.3_reduced_TCA_genome.fa   \
  --annotationA testset/NC_035109.1_reduced_AAE_genome.gff   \
  --annotationB testset/NC_007416.3_reduced_TCA_genome.gff   \
  --clip testset/SRR5163632_aae_clip_reduced.fastq,testset/SRR5163633_aae_clip_reduced.fastq,testset/SRR5163634_aae_clip_reduced.fastq   \
  --clip testset/SRR5163635_aae_clip_reduced.fastq,testset/SRR5163636_aae_clip_reduced.fastq,testset/SRR5163637_aae_clip_reduced.fastq --adapterclip GTGTCAGTCACTTCCAGCGG  \
  --overwrite \
  --smallrnaseq a=testset/tca_smallRNAseq_rna_contaminated.fastq \
  --adaptersmallrnaseq3=TGGAATTCTCGGGTGCCAAGG \
  --adaptersmallrnaseq5 GTTCAGAGTTCTACAGTCCGACGATC \
  --filterncrnas testset/TCA_all_ncRNA_but_miR.fa \
  --speciesB tca 2>&1 | tee out.log
```


## Usage

### Input data
  - minimal workflow
    - speciesA genome
    - speciesA GFF
    - speicesA AGO-CLIP-sequencing library/libraries
    - speciesB genome
    - speciesB GFF
    - speciesB microRNA set (mature)
  - full workflow (in addition to the minimal workflow)
    - speciesB non-codingRNA set (without miRNAs)
    - speciesB microRNA set (precursor)
    - speciesB smallRNA-sequencing library/libraries
    
### PARAMETERS

- `--version|-V`

    version of this pipeline

- `--help|-h`

    prints a helpful help message

- `--genomeA` and `--genomeB`

    Genome of the species with the CLIP data (species A, `--genomeA`) and
    the genome of the species where we want to predict the miRNA targets
    (species B, `--genomeB`)

- `--gffA` and `--gffB`

    Genome feature file (GFF) of the species with the CLIP data (species
    A, `--gffA`) and the GFF of the species where we want to predict the
    miRNA targets (species B, `--gffB`)

- `--clip`

    Comma-separated CLIP-seq .fastq files in Format

        --clip con1_rep1_clip.fq,con1_rep2_clip.fq,con2_clip.fq
        # OR
        --clip con1_rep1_clip.fq --clip con1_rep2_clip.fq --clip con2_clip.fq

- `--adapterclip`

    Sequencing-adapter of CLIP reads

- `--smallrnaseq`

    Comma-separated smallRNA-seq FASTQ files, initialized with
    'condition=' in Format

        --smallrnaseq con1=A.fastq,B.fastq --smallrnaseq con2=C.fq
        # OR
        --smallrnaseq con1=A.fastq --smallrnaseq con1=B.fastq --smallrnaseq con2=C.fq

- `--adaptersmallrnaseq5` and `--adaptersmallrnaseq3`

    5' adapter of smallRNA-seq reads (`--adaptersmallrnaseq5`) and for 3' end (`--adaptersmallrnaseq3`)

- `--filterncrnas`

    Multi-fasta file of ncRNAs to filter smallRNA-seq reads. Those must
    not contain miRNAs.

- `--threads`

    Number of threads to be used

- `--overwrite`

    set this parameter to overwrite existing files

- `--testrun`

    sets this pipeline to testmode (accounting for small testset in
    piranha). This option should not be used in real analysis!

- `--out`

    output folder

- `--mirna`

    miRNA set, if set, mining is disabled and this set is used for prediction

- `--speciesBtag`

    Three letter code of species where we want to predict the miRNA
    targets (species B, `--speciesBtag`).

- `--mirbasedir`

    The folder specified by `--mirbasedir` is searched for the files
    `organisms.txt.gz`, `mature.fa.gz`, and `hairpin.fa.gz`. If the
    files are not exist, they will be downloaded.

- `--tempdir`

    The folder specified by `--tempdir` is used for temporary files. The
    default value is `tmp/` inside the output folder specified by the
    `--out` parameter.

- `--piranahbinsize`

    Sets the `Piranah` bin size and has a default value of `20`.

### OUTPUT

- pseudo mirBASE dat file: `final_mirbase_pseudofile.dat`

    A pseudo mirBASE dat file containing all precursor sequences with their named mature sequences and their coordinates. It only contain the fields:

    - `ID`
    - `FH` and `FT`
    - `SQ`

- mature miRNA set: `mature_combined_mirbase_novel.fa`

    mature microRNA set, containing novels and miRBase-completed (if mined), together with the known miRNAs from miRBase

- precursor miRNA set: `hairpin_combined_mirbase_novel.fa`

    precursor microRNA set, containing novels (if mined), together with the known miRNAs from miRBase

- mature miRNA expression per condition: `miRNA_expression.csv`

    Semicolon-separated file containing:

    - 1. `rpm`
    - 2. `condition`
    - 3. `miRNA`

- orthologous prediction file: `miRNA_orthologs.csv`

    tab-separated file containing:

    - 1. `query_id`
    - 2. `subject_id`
    - 3. `identity`
    - 4. `alignment length`
    - 5. `number mismatches`
    - 6. `number gap openings`
    - 7. `start position inside query`
    - 8. `end position inside query`
    - 9. `start position inside subject`
    - 10. `end position inside subject`
    - 11. `evalue`
    - 12. `bitscore`
    - 13. `aligned query sequence`
    - 14. `aligned subject sequence`
    - 15. `length query sequence`
    - 16. `length subject sequence`
    - 17. `coverage for query sequence`
    - 18. `coverage for subject sequence`

- miRDeep2 mining result in HTML/CSV `mirdeep_output.html/csv`

    the standard output HTML/CSV file of miRDeep2

- ISOMIR prediction files: `isomir_output_CONDITION.csv`

    semincolon delimited file containing:

    - 1. `mirna`
    - 2. `substitutions`
    - 3. `added nucleotids on 3' end`
    - 4. `nucleotides at 5' end different from the annonated sequence`
    - 5. `nucleotides at 3' end different from the annonated sequence`
    - 6. `sequence`
    - 7. `rpm`
    - 8. `condition`

- genomics location of miRNAs: `miRNA_genomic_position.csv`

    tab delimited file containing:

    - 1. `miRNA`
    - 2. `genomic contig`
    - 3. `identify`
    - 4. `length`
    - 5. `miRNA-length`
    - 6. `number mismatches`
    - 7. `number gapopens`
    - 8. `miRNA-start`
    - 9. `miRNA-stop`
    - 10. `genomic-start`
    - 11. `genomic-stop`
    - 12. `evalue`
    - 13. `bitscore`

- all library support-level target predictions: `*_miranda_output.txt`

    miranda output, reduced to the lines, starting with > only

- all library support-level CLIP transfer .bed files: `*transfered_merged.bed`

    bed-file of the transferred CLIP-regions in speciesB transcriptome

### Example
#### Testset
Feel free to test the pipeline with our [microPIECE-testset](https://github.com/microPIECE-team/microPIECE-testset) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1188471.svg)](https://doi.org/10.5281/zenodo.1188471):

`git clone git@github.com:microPIECE-team/microPIECE-testset.git`

#### Alternative
  - **minimal workflow**
    - speciesA genome [AAE genome : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz](https://tinyurl.com/yagl5mlo)
    - speciesA GFF [AAE GFF : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.gff.gz](https://tinyurl.com/ybckl5pp)
    - speicesA AGO-CLIP-sequencing library/libraries [AGO-CLIP of AAE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93345)
    - speciesB genome [TCA genome : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.fna.gz](https://tinyurl.com/y7844w3t)
    - speciesB GFF [TCA GFF : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.gff.gz](https://tinyurl.com/ybj35n7j)
    - speciesB microRNA set (mature) [TCA mature miRNAs](http://mirbase.org/cgi-bin/mirna_summary.pl?org=tca)
    
  - **full workflow** (in addition to the minimal workflow)
    - speciesB microRNA set (precursor) [TCA stem-loop miRNAs](http://mirbase.org/cgi-bin/mirna_summary.pl?org=tca)
    - speciesB smallRNA-sequencing library/libraries [TCA smallRNA-sequencing data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63770)
    - speciesB non-codingRNA set (without miRNAs - to filter smRNA-seq data) (*OPTIONAL*)

## CAVEATS
Complete list of open issues is available on [Github-Issues](https://github.com/microPIECE-team/microPIECE/issues).

Please report any new issues ad [new Github-Issue](https://github.com/microPIECE-team/microPIECE/issues/new).

## Changelog
- scheduled for next release

    No features planed

- [v1.4.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.4.0) (2018-03-31)

    Copying pseudo mirBASE dat file `final_mirbase_pseudofile.dat` into output folder (Fixes [#131](https://github.com/microPIECE-team/microPIECE/issues/131))

    Corrected `RNA::HairpinFigure` output (Fixes [#137](https://github.com/microPIECE-team/microPIECE/issues/137))

    Fix the requirement of an accession inside mirBASE dat file (Fixes [#134](https://github.com/microPIECE-team/microPIECE/issues/134))

    Avoiding error message while copying the out file for genomic location into base folder (Fixes [#117](https://github.com/microPIECE-team/microPIECE/issues/117)) 

- [v1.3.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.3.0) (2018-03-29)

    Creating all structures on the fly using pseudo-mirBASE-dat as input.

    Using `miRNA.dat` from mirBASE as source for mature/precursor sequence and relationship (Fixes [#127](https://github.com/microPIECE-team/microPIECE/issues/127))

    Fix of division-by-zero bug for empty mapping files (Fixes [#118](https://github.com/microPIECE-team/microPIECE/issues/118))

    Fix of typo in `--piranhabinsize` option (Fixes [#116](https://github.com/microPIECE-team/microPIECE/issues/116))

- [v1.2.3](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.3) (2018-03-26)

    Fix transformation of precursor sequences based on mirbase #22 precursor sequences with a single mature.
    (Fixes L<#109|https://github.com/microPIECE-team/microPIECE/issues/109>)

- [v1.2.2](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.2) (2018-03-23)

    Improved collision detection for newly identified miRNAs avoiding crashed caused by genomic copies.
    (Fixes [#105](https://github.com/microPIECE-team/microPIECE/issues/105))

- [v1.2.1](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.1) (2018-03-23)

    Enables stable numbering for newly identified miRNAs based on their precursor and mature sequences
    (Fixes [#101](https://github.com/microPIECE-team/microPIECE/issues/101))

- [v1.2.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.2.0) (2018-03-22)

    We are using miraligner which requires a java version 1.7, but 1.8 was
    installed by default. This was fixed by switching to v1.4 of the
    docker base image. Additionally, miraligner requires fix filenames for
    its databases. Therefore, the version v1.2.0 solved miraligner related
    bugs and reenables the isomir detection.  (Fixes
    [#97](https://github.com/microPIECE-team/microPIECE/issues/97) and
    [#98](https://github.com/microPIECE-team/microPIECE/issues/98))

- [v1.1.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.1.0) (2018-03-12)

    Add isomir detection and copy the final genomic location file to the
    output filter (Fixes
    [#34](https://github.com/microPIECE-team/microPIECE/issues/34))

- [v1.0.7](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.7) (2018-03-08)

    Piranha was lacking of a bin\_size parameter. Added parameter `--piranahbinsize` with a default value of `20`
    (Fixes [#66](https://github.com/microPIECE-team/microPIECE/issues/80))

- [v1.0.6](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.6) (2018-03-08)

    Added parameter `--mirbasedir` and `--tempdir` to support local
    mirbase files and relocation of directory for temporary files (Fixes
    [#66](https://github.com/microPIECE-team/microPIECE/issues/66),
    [#73](https://github.com/microPIECE-team/microPIECE/issues/73), and
    [#76](https://github.com/microPIECE-team/microPIECE/issues/76))

- [v1.0.5](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.5) (2018-03-07)

    Update of documentation and correct spelling of `--mirna` parameter

- [v1.0.4](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.4) (2018-03-07)

    Fixes complete mature in final output (Fixes [#69](https://github.com/microPIECE-team/microPIECE/issues/69))

- [v1.0.3](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.3) (2018-03-06)

    Add tests for perl scripts in script folder which ensure the correct handling of BED stop coordinates (Fixes [#65](https://github.com/microPIECE-team/microPIECE/issues/65))

- [v1.0.2](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.2) (2018-03-05)

    Fixes the incorrect sorting of BED files, result was correct, but sorting was performed in the wrong order. (Fixes [#63](https://github.com/microPIECE-team/microPIECE/issues/63))

- [v1.0.1](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.1) (2018-03-05)

    Fix an error conserning BED file handling of start and stop coordinates. (Fixes [#59](https://github.com/microPIECE-team/microPIECE/issues/59))

- [v1.0.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v1.0.0) (2018-03-05)

    <div>
            is archived as <a href="https://doi.org/10.5281/zenodo.1188484"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1188484.svg" alt="DOI" /></a>
            and submitted to <a href="http://joss.theoj.org">The Journal of Open Source Software</a>.
    </div>

- [v0.9.0](https://github.com/microPIECE-team/microPIECE/releases/tag/v0.9.0) (2018-03-05)

    <div>
            first version archived at Zenodo with the <a href="https://doi.org/10.5281/zenodo.1188481"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1188481.svg" alt="DOI" /></a>
    </div>

## License
This program is released under GPLv2. For further license information, see LICENSE.md shipped with this program.
Copyright(c)2018 Daniel Amsel and Frank Förster (employees of Fraunhofer Institute for Molecular Biology and Applied Ecology IME) All rights reserved.

# AUTHORS
- Daniel Amsel &lt;daniel.amsel@ime.fraunhofer.de>
- Frank Förster &lt;frank.foerster@ime.fraunhofer.de>

# SEE ALSO
[Project source code on Github](https://github.com/microPIECE-team/microPIECE)
[Docker image on DockerHub](https://hub.docker.com/r/micropiece/micropiece/)
[Travis continuous integration page](https://travis-ci.org/microPIECE-team/microPIECE)
[Test coverage reports](https://coveralls.io/github/microPIECE-team/microPIECE)
