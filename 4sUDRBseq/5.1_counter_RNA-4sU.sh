#!/bin/bash

# This script counts the number of aligned reads per sample

SORTEDBAM="/path/4sUDRB/Alignments/*.sam"

for FILE in $SORTEDBAM
do	
	NUMBER=$(cat $FILE | grep -vP "^@" | wc -l)
	NAME="$(basename $FILE)"
	NAME=$(echo $NAME | sed 's/\.sam//')
	printf "$NAME\t$NUMBER\n" >> READ_NUMBER.txt
	echo $NAME
done
