#!/bin/bash

# This script sorts alignments ordered positionally based upon their alignment coordinates on each chromosome

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	BAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/1_SAMtoBAM/MG9-"$NUM".bam"
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/MG9-"$NUM".sorted.bam"

	samtools sort $BAM -o $SALIDA
done
wait
