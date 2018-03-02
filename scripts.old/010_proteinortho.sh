#!/bin/bash
command -v makeblastdb >/dev/null 2>&1 || { echo "I require makeblastdb but it's not installed. Aborting." >&2; exit 1; }
command -v proteinortho5.pl >/dev/null 2>&1 || { echo "I require proteinortho5.pl but it's not installed. Aborting." >&2; exit 1; }

#create folder for this step
mkdir ../010/
#unzip data
zcat ../data/TCA/GCF_000002335.3_Tcas5.2_protein.faa.gz >../010/GCF_000002335.3_Tcas5.2_protein.faa
zcat ../data/AAE/GCF_000004015.4_AaegL3_protein.faa.gz  >../010/GCF_000004015.4_AaegL3_protein.faa
# makeblastdb
makeblastdb -in ../010/GCF_000002335.3_Tcas5.2_protein.faa -dbtype prot
makeblastdb -in ../010/GCF_000004015.4_AaegL3_protein.faa -dbtype prot
# run proteinortho
/opt/proteinortho_v5.16/proteinortho5.pl  -clean --project=../010/TCA_vs_AAE ../010/GCF_000002335.3_Tcas5.2_protein.faa ../010/GCF_000004015.4_AaegL3_protein.faa
rm ../010/GCF_000002335.3_Tcas5.2_protein.faa
rm ../010/GCF_000004015.4_AaegL3_protein.faa
