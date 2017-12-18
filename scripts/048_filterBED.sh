#!/bin/bash
mkdir -p ../040/

for i in 1 2 3 4 5 6;
do
    echo "Working on setting $i of 6"
    ./049_bed2signal.pl ../040/clip_merged.bed "$i" > ../040/clip_merged_"$i"of6BEDfilter.bed
done


