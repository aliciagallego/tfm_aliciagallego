#!/bin/bash

BAM="/path/cheRNA/Pull/Alignment/*.bam"
OUT="/path/cheRNA/Pull"

for FILE in $BAM
do	
	NUMBER=$(samtools view $FILE | wc -l)
	NAME="$(basename $FILE)"
	NAME=$(echo $NAME | sed 's/\.bam//')
	printf "$NAME\t$NUMBER\n" >> $OUT/READ_NUMBER_pull.txt
done
