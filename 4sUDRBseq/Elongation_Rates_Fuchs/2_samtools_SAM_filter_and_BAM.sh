#!/bin/bash

# This script filters SAM alignments Skip alignments with MAPQ smaller than 1 and converts SAM format into BAM

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	SAM="/media/cc/B/Josemi/TTseq_Feb2022/Alignments/MG10_RNA_DRB-4sU/MG9-"$NUM".sam"
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/6_SAM_filtered_q1/MG9-"$NUM"_filtered.bam"
	
	# despues cambie a mano los outputs por .sam (creo que lo que se generan son SAM no BAM)

	samtools view -S -q 1 -b $SAM > $SALIDA
done
wait



