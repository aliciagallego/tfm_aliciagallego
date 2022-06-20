#!/bin/bash

# This script converts BAM files to SAM 

for file in /media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/Pull/*bam
do
	
	NAMEE=$(basename $file)
	NAME=${NAMEE%.*}
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/Pull/Pull_sam/"$NAME".sam"
    	echo $file
    	samtools view -h $file > $SALIDA
done
