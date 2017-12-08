# run gsnap to map the reads of CLIP sequencing to AAE genome
command -v gsnap >/dev/null 2>&1 || { echo "I require gsnap but it's not installed. Aborting." >&2; exit 1; }
command -v gmap_build >/dev/null 2>&1 || { echo "I require gsnap but it's not installed. Aborting." >&2; exit 1; }

gmap_build -D gsnap_db -k 15 -d Aedes_aegypti GCF_000004015.4_AaegL3_genomic.fna

