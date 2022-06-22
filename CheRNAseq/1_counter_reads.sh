#!/bin/bash

BAM="/path/cheRNA/Alignments/*.bam"
OUT="/path/cheRNA/"

for FILE in $BAM
do	
	NUMBER=$(samtools view $FILE | wc -l)
	NAME="$(basename $FILE)"
	NAME=$(echo $NAME | sed 's/\.bam//')
	printf "$NAME\t$NUMBER\n" >> $OUT/READ_NUMBER.txt
done
