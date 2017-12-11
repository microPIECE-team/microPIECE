#!/bin/bash
command -v ./4_map_clip_gff_needle.pl >/dev/null 2>&1 || { echo "I require 4_map_clip_gff_needle.pl but it's not installed. Aborting." >&2; exit 1; }
mkdir ../080/
# filter GFF files for unqiue longest transcripts
zcat ../data/AAE/GCF_000004015.4_AaegL3_genomic.gff.gz > ../080/GCF_000004015.4_AaegL3_genomic.gff
# AAE
perl -aF'/\t/' -ne 'next unless ($F[2] =~ /mRNA|exon|CDS/); $F[8] =~ /Dbxref=VectorBase:([^-]+)-[RP]([^,]+).+Genbank:([^,]+)/; if ($F[2] eq "mRNA") { $h{$1}{$2} = {}; } elsif ($F[2] eq "exon" && exists $h{$1}{$2} ) { $h{$1}{$2}{exon}{name} = $3; $h{$1}{$2}{exon}{len} += abs($F[3]-$F[4]); } elsif (exists $h{$1}{$2}) { $h{$1}{$2}{cds}{name} = $3; } else { warn "For line $_ no mRNA was found\n"; } END{ foreach my $vec(sort keys %h) { my ($longest) = sort {$h{$vec}{$b}{exon}{len} <=> $h{$vec}{$a}{exon}{len}} (keys %{$h{$vec}}); next unless (exists $h{$vec}{$longest}{exon} && exists $h{$vec}{$longest}{cds}); print join("\t", sort ($h{$vec}{$longest}{cds}{name}, $h{$vec}{$longest}{exon}{name})), "\n"; }}' -MData::Dumper ../080/GCF_000004015.4_AaegL3_genomic.gff > ../080/GCF_000004015.4_AaegL3_genomic_XM_XP_unique.csv 2> ../080/GCF_000004015.4_AaegL3_genomic_XM_XP_unique.err
# TCA

zcat ../data/TCA/GCF_000002335.3_Tcas5.2_rna.fna.gz > ../080/GCF_000002335.3_Tcas5.2_rna.fna

# transfer with needle
./4_map_clip_gff_needle.pl ../080/GCF_000004015.4_AaegL3_genomic_XM_XP_unique.csv ../010/TCA_vs_AAE.proteinortho ../070/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC.fa ../080/GCF_000002335.3_Tcas5.2_genomic_XM_XP_unique.csv ../080/GCF_000002335.3_Tcas5.2_rna.fna ../080/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle.csv > ../080/SRR5163632_SRR5163633_SRR5163634_SRR5163635_SRR5163636_SRR5163637_trim_gsnap_piranha_sort_merge_mapGFF_minLen0_intersect_min22_max50_cat_sort_merge_UC_mapNeedle.aln

