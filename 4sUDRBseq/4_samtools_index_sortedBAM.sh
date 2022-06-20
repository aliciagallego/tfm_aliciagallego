#!/bin/bash

# This script indexes sorted BAM files to quickly extract alignments overlapping particular genomic regions and view alignments in a genome browser such as IGV
# -b creates a BAI index (this is currently the default when no format options are used)

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	SORTEDBAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/MG9-"$NUM".sorted.bam"

	samtools index $SORTEDBAM
done
wait
