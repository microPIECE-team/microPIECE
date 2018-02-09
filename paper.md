---
title: 'miRNA CLIP pipe'
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
- name: Frank Förster
  orcid: 0000-0003-4166-5423
  affiliation: 1
- name: Andreas Vilcinskas
  orcid:
  affiliation: "1, 2"
- name: André Billion
  orcid:
  affiliation: 1
affiliations:
- name: Fraunhofer Institute for Molecular Biology and Applied Ecology, Department of Bioresources, Winchester Str. 2, 35394 Giessen, Germany
  index: 1
- name: Institute for Insect Biotechnology, Heinrich-Buff-Ring 26-32, 35392 Giessen, Germany
  index: 2
date: 06 February 2018
bibliography: paper.bib
---

# Summary
MicroRNAs are post-transcriptional fine-regulators. With a length of around 21 nucleotides, they form a RNA-induced silencing complex (RISC) complex with a protein of the Argonaute family. This complex then binds to the mRNA UTR and CDS regions and in general promotes degradation or translational inhibition. It is now interesting to know the microRNA-mRNA pairs in order to infer dysregulating effects on the organism. In order to assign a microRNA to a mRNA target, various tools with different technical approaches were developed. They are mostly based on the assumption that the first eight nucleotides of the microRNA (seed region) determines the binding region on the mRNA. Some approaches also include supporting bindings in the rear part of the microRNA, others include secondary structures of the mRNA or binding energies of the mRNA-miRNA complex. But they all suffer from the statistical problem that such short sequences, like a microRNAs or even the seed regions, can occur simply by chance. This results in a large amount of false positive predictions. Therefore, wet lab validation methods like luciferase reporter assays were adapted. The disadvantage here is that those methods are not applicable for high-throughput analysis. Another, more scalable method is cross-linking immunoprecipitation-high-throughput sequencing (CLIP).  Here, binding regions of the RISC show a specific signal that can be used to shrink the search space of miRNA target predictions. The limitation here is the difficult technical progress in the laboratory. This is the reason why there are only a few datasets available for human, mouse, worm and mosquito. 


- A summary describing the high-level functionality and purpose of the software
for a diverse, non-specialist audience
- A clear statement of need that illustrates the purpose of the software
- A list of key references including a link to the software archive
- Mentions (if applicable) of any ongoing research projects using the software
or recent scholarly publications enabled by it

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

This is an example citation [@figshare_archive].

Figures can be included like this: ![Fidgit deposited in figshare.](figshare_article.png)

# References
