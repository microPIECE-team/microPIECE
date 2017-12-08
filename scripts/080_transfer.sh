#!/bin/bash
command -v ./4_map_clip_gff_needle.pl >/dev/null 2>&1 || { echo "I require 4_map_clip_gff_needle.pl but it's not installed. Aborting." >&2; exit 1; }

# filter GFF files for unqiue longest transcripts
# AAE
perl -aF'/\t/' -ne 'next unless ($F[2] =~ /mRNA|exon|CDS/); $F[8] =~ /Dbxref=VectorBase:([^-]+)-[RP]([^,]+).+Genbank:([^,]+)/; if ($F[2] eq "mRNA") { $h{$1}{$2} = {}; } elsif ($F[2] eq "exon" && exists $h{$1}{$2} ) { $h{$1}{$2}{exon}{name} = $3; $h{$1}{$2}{exon}{len} += abs($F[3]-$F[4]); } elsif (exists $h{$1}{$2}) { $h{$1}{$2}{cds}{name} = $3; } else { warn "For line $_ no mRNA was found\n"; } END{ foreach my $vec(sort keys %h) { my ($longest) = sort {$h{$vec}{$b}{exon}{len} <=> $h{$vec}{$a}{exon}{len}} (keys %{$h{$vec}}); next unless (exists $h{$vec}{$longest}{exon} && exists $h{$vec}{$longest}{cds}); print join("\t", sort ($h{$vec}{$longest}{cds}{name}, $h{$vec}{$longest}{exon}{name})), "\n"; }}' -MData::Dumper GCF_000004015.4_AaegL3_genomic.gff > GCF_000004015.4_AaegL3_genomic_XM_XP_unique.csv

# TCA

