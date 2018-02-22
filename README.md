# microPIECE

The `microPIECE` (microRNA pipeline enhanced by CLIP experiments) takes the AGO-CLIP data from a speciesA and transfers it to a speciesB. Given a set of miRNAs from speciesB it then predicts their targets on the transfered CLIP regions.

For the minimal workflow it needs a genome file, as well as its annotation file in GFF format for speciesA and speciesB. For speciesA at least one AGO-CLIP dataset is needed and speciesB needs a set of miRNAs for the target prediction. For the full workflow, a set of smallRNA-sequencing data is additionally needed and a set of non-coding RNAs can be provided as filter. The pipeline uses the smallRNA data for the mining of novel microRNAs and the completion of the given miRNA dataset, if needed. It further performs expression calculation, isoform detection, genomic loci identification and orthology determination.

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
`git clone git@github.com:DanielAmsel/microPIECE.git`

## Docker
We also provide `microPIECE` as DOCKER image. We tested the image on Ubuntu, Debian and MacOS. For the latter one, the `Piranha` command `make test` fails during the build, but when entering the container, the test succeds. Therefore, we temporarily excluded this statement.

## Usage

## Input data
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

## Example

  - **minimal workflow**
    - speciesA genome [AAE genome](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz)
    - speciesA GFF [AAE GFF](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.gff.gz)
    - speicesA AGO-CLIP-sequencing library/libraries [AGO-CLIP of AAE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93345)
    - speciesB genome [TCA genome](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.fna.gz)
    - speciesB GFF [TCA GFF] (ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/335/GCF_000002335.3_Tcas5.2/GCF_000002335.3_Tcas5.2_genomic.gff.gz)
    - speciesB microRNA set (mature) [TCA mature miRNAs](http://mirbase.org/cgi-bin/mirna_summary.pl?org=tca)
    
  - **full workflow** (in addition to the minimal workflow)
    - speciesB microRNA set (precursor) [TCA stem-loop miRNAs](http://mirbase.org/cgi-bin/mirna_summary.pl?org=tca)
    - speciesB smallRNA-sequencing library/libraries [TCA smallRNA-sequencing data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63770)
    - speciesB non-codingRNA set (without miRNAs - to filter smRNA-seq data) (OPTIONAL)

    
## Changelog
Version 1.0.0 is archived as *DOI* and submitted to [The Journal of Open Source Software](http://joss.theoj.org/).
## License
This program is released under GPLv2. For further license information, see LICENSE.md shipped with this program.
Copyright(c)2018 Daniel Amsel and Frank Förster (employees of Fraunhofer Institute for Molecular Biology and Applied Ecology IME) All rights reserved.
