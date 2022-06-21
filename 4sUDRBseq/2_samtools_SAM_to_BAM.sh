#!/bin/bash

# This script converts SAM format into BAM

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	SAM="/path/4sUDRB/MG9-"$NUM".sam"
	SALIDA="/path/4sUDRB/MG9-"$NUM".bam"

	samtools view -S -b $SAM > $SALIDA
done
wait
