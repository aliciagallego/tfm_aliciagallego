#!/bin/bash

# This script counts the number of aligned reads

SORTEDBAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/*.sorted.bam"

for FILE in $SORTEDBAM
do
	NUMBER=$(samtools idxstats $FILE | awk 'BEGIN {sum=0}; {sum+=$3}; END {print sum}')
	NAME="$(basename $FILE)"
	NAME=$(echo $NAME | sed 's/\.sorted.bam//')
	printf "$NAME\t$NUMBER\n" >> READ_NUMBER.txt
done
