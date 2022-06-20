#!/bin/bash

# This script counts the number of aligned reads per sample

SORTEDBAM="/media/cc/B/Josemi/TTseq_Feb2022/Alignments/MG10_RNA_DRB-4sU/*.sam"

for FILE in $SORTEDBAM
do	
	NUMBER=$(cat $FILE | grep -vP "^@" | wc -l)
	NAME="$(basename $FILE)"
	NAME=$(echo $NAME | sed 's/\.sam//')
	printf "$NAME\t$NUMBER\n" >> READ_NUMBER.txt
	echo $NAME
done
