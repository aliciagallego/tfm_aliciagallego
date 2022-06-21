#!/bin/bash

# This script converts BAM files to SAM 

for file in /path/4suDRB/Pull/*bam
do
	NAMEE=$(basename $file)
	NAME=${NAMEE%.*}
	SALIDA="/path/4suDRB/Pull/"$NAME".sam"
    	echo $file
    	samtools view -h $file > $SALIDA
done
