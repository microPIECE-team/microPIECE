#!/usr/bin/bash
#create folder for this step
mkdir ../010/
#unzip data
zcat ../data/GCF_000002335.3_Tcas5.2_protein.faa.gz >../010/GCF_000002335.3_Tcas5.2_protein.faa
zcat ../data/GCF_000004015.4_AaegL3_protein.faa.gz  >../010/GCF_000004015.4_AaegL3_protein.faa
# makeblastdb
makeblastdb ../010/GCF_000002335.3_Tcas5.2_protein.faa
makeblastdb ../010/GCF_000004015.4_AaegL3_protein.faa
# run proteinortho
proteinortho5.pl ../010/GCF_000002335.3_Tcas5.2_protein.faa ../010/GCF_000004015.4_AaegL3_protein.faa
