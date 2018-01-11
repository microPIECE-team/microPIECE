bwa index 061_mature.fa 
for i in *.fq; do bwa aln -n 1 -o 0 -e 0 -k 1 -t 10 -f ${i}.sai 061_mature.fa $i; done
for i in *.fq; do bwa samse -f ${i}.sam 061_mature.fa ${i}.sai $i; done
for i in *.fq; do samtools view -F 4 ${i}.sam -o ${i}_aln.sam; done
for i in *.fq; do ./xa2multi.pl ${i}_aln.sam > ${i}_aln_xa2multi.sam; done
