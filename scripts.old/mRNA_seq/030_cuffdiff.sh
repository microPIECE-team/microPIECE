#!/bin/bash
mkdir -p data/030_cuffdiff/ 
cuffdiff -o data/030_cuffdiff/ -L oocyte,L5,LLI db/GCF_000002335.3_Tcas5.2_genomic.gff data/020_mapped_RNA/oocyte/accepted_hits.bam data/020_mapped_RNA/L5/accepted_hits.bam data/020_mapped_RNA/LLI/accepted_hits.bam
