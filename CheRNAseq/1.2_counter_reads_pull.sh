#!/bin/bash

BAM="/media/cc/A/Josemi/NGS/cheRNA/BAMs/Pulls/*.bam"
OUT="/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output"

for FILE in $BAM
do	
	NUMBER=$(samtools view $FILE | wc -l)
	NAME="$(basename $FILE)"
	NAME=$(echo $NAME | sed 's/\.bam//')
	printf "$NAME\t$NUMBER\n" >> $OUT/READ_NUMBER_pull.txt
done
