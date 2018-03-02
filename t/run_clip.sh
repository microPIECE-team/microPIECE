#!/bin/bash
cd /opt/microPIECE

perl -MDevel::Cover=-db,coverage,statement,branch,condition,path,subroutine,time /opt/microPIECE/microPIECE.pl \
   --genomeA /tmp/microPIECE-testset/NC_035109.1_reduced_AAE_genome.fa \
   --genomeB /tmp/microPIECE-testset/NC_007416.3_reduced_TCA_genome.fa \
   --annotationA /tmp/microPIECE-testset/NC_035109.1_reduced_AAE_genome.gff \
   --annotationB /tmp/microPIECE-testset/NC_007416.3_reduced_TCA_genome.gff \
   --clip /tmp/microPIECE-testset/SRR5163632_aae_clip_reduced.fastq,/tmp/microPIECE-testset/SRR5163633_aae_clip_reduced.fastq,/tmp/microPIECE-testset/SRR5163634_aae_clip_reduced.fastq \
   --clip /tmp/microPIECE-testset/SRR5163635_aae_clip_reduced.fastq,/tmp/microPIECE-testset/SRR5163636_aae_clip_reduced.fastq,/tmp/microPIECE-testset/SRR5163637_aae_clip_reduced.fastq \
   --adapterclip GTGTCAGTCACTTCCAGCGG \
   --mirna /tmp/microPIECE-testset/tca_known_mirs.fa \
   --out clip
