#!/bin/bash

# This script converts SAM format into BAM

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	SAM="/media/cc/B/Josemi/TTseq_Feb2022/Alignments/MG10_RNA_DRB-4sU/MG9-"$NUM".sam"
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/1_SAMtoBAM/MG9-"$NUM".bam"

	samtools view -S -b $SAM > $SALIDA
done
wait
