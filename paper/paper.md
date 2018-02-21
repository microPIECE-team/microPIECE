---
title: 'microRNA CLIP pipe'
tags:
- CLIP-seq
- microRNA
- target prediction
- homology
- insectomics
authors:
- name: Daniel Amsel
  orcid: 0000-0002-0512-9802
  affiliation: 1
- name: André Billion
  orcid:
  affiliation: 1
- name: Andreas Vilcinskas
  orcid: 0000-0001-8276-4968
  affiliation: "1, 2"
- name: Frank Förster
  orcid: 0000-0003-4166-5423
  affiliation: 1
affiliations:
- name: Fraunhofer Institute for Molecular Biology and Applied Ecology, Department of Bioresources, Winchester Str. 2, 35394 Giessen, Germany
  index: 1
- name: Institute for Insect Biotechnology, Heinrich-Buff-Ring 26-32, 35392 Giessen, Germany
  index: 2
date: 21 February 2018
bibliography: paper.bib
---

# Summary
All microRNAs are assumed to be post-transcriptional fine-regulators. With a length of around 21 nucleotides, they form a RNA-induced silencing complex (RISC) complex with a protein of the Argonaute family. This complex then binds to the messengerRNA untranslated regions and coding sequence regions and in general promotes degradation or translational inhibition. It is now important to know the microRNA-mRNA pairs in order to infer dysregulating effects on the organism. In order to assign a microRNA to a mRNA target, various tools with different technical approaches were developed. They are mostly based on the assumption that the first eight nucleotides of the microRNA (seed region) determine the binding region on the mRNA. Some approaches also include supporting bindings in the rear part of the microRNA, others take secondary structures of the mRNA or binding energies of the mRNA-miRNA complex into account. Nevertheless, they all suffer from the statistical problem that such short target regions, often occur simply by chance in transcript sequences. This results in a huge amount of false positive predictions. A target prediction of all 430 *Tribolium castaneum* mature microRNAs from <miRBase.org> against all 18.534 protein coding cDNA sequences from <Ensembl.org> results in 2.203.593 possible microRNA-target interactions, predicted by the commonly used tool `miranda` [@betel2008microrna] with standard parameters. To increase the credibility, wet lab validation methods like luciferase reporter assays are required. The disadvantage here is that this workflow is not applicable for high-throughput analysis, as it can only treat small subsets of sequence combinations. Another, more scalable method is cross-linking immunoprecipitation-high-throughput sequencing (CLIP-seq). Here, binding regions of the RISC show a specific signal in the sequencing reads that can be used to shrink the search space of miRNA target predictions, when mapping them to the transcriptome. The limitation here is the difficult technical treatment in the laboratory. This is the reason why there are only a few datasets available for human, mouse, worm and mosquito. It would now be useful, if we could simply transfer the information of a binding region, already identified by CLIP-seq, to another species. This is what our microRNA CLIP pipeline is about. 

The pipeline (Figure 1) starts with the *speciesA* CLIP-seq library trimming, using cutadapt [@martin2011cutadapt]. Trimmed reads are then mapped to the genome with gsnap [@wu2005gmap] and the results are evaluated by Piranha [@uren2012site]. Then the libraries are merged into the `BED` file format. This file includes a column that displays how many libraries support each genomic position. Next, a file for each library-support-level is created, so that the user can decide how many CLIP libraries are necessary to account this region as binding region. Now for each library-support-level, an assignment of each sequence to the transcriptome is performed. Outgoing from the transcript, the corresponding protein is used to discover the orthologous protein in the *speciesB* by proteinortho [@lechner2011proteinortho]. This information is used as criteria to align the CLIP region from *speciesA* to the orthologous transcript in *speciesB* with needle [@rice2000emboss].

The miRNA analysis part can either treat a predefined set of miRNAs or can additionally take smallRNA sequencing libraries into account. In the latter case, it uses cutadapt [@martin2011cutadapt] to trim the adapter sequences from the small RNA sequencing libraries from *speciesB*. The trimmed libraries are then filtered for ncRNAs using bwa [@li2009fast].The resulting files are merged into a pooled set and used for mining of novel microRNAs with miRDeep2 [@friedlander2011mirdeep2]. The pipeline then parses the result file and tries to add missing entries from the <miRBase.org> database [@kozomara2013mirbase], e.g. if only one arm was previously annotated and the mining discoveres the exact position of the particular arm. Afterwards the expression of each miRNA, including the novel miRNAs, is calculated as ReadsPerMillion (RPM), outgoing from the non-pooled trimmed and filtered libraries. The pipeline also accounts for miRNA isoforms (isomiRs), by removing all trimmed reads, containing undetermined nucleotides. It then passes the reads to miraligner from the seqbuster package [@pantano2009seqbuster]. We then perform a target prediction with `miranda` [@betel2008microrna] on the previously transfered orthologous CLIP regions. Orthologous miRNAs in other species were determined by a `BLASTN` [@altschul1990basic] search against all metazoan miRNAs from <miRBase.org>. Finally, the genomic regions for the miRNAs were also identified by a `BLASTN` search against the genome of *speciesB*.

![Scheme of the pipeline](miRNA_CLIP_pipe.png)

We further used samtools [@li2009sequence] and bedtools [@quinlan2010bedtools] for file conversion.

# References
