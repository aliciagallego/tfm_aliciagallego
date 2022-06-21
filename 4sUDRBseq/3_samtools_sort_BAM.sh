#!/bin/bash

# This script sorts alignments ordered positionally based upon their alignment coordinates on each chromosome

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	BAM="/path/4sUDRB/MG9-"$NUM".bam"
	SALIDA="/path/4sUDRB/MG9-"$NUM".sorted.bam"

	samtools sort $BAM -o $SALIDA
done
wait
