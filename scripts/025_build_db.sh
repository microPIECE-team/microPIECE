#!/bin/bash
# creates a GSNAP database

command -v gmap_build >/dev/null 2>&1 || { echo "I require gmap_build but it's not installed. Aborting." >&2; exit 1; }



# creation of temporary folder for this step
mkdir ../030/

# build genome index for gsnap
gmap_build -g -D ../030/gsnap_db -k 15 -d Aedes_aegypti ../data/AAE/GCF_000004015.4_AaegL3_genomic.fna.gz

