#!/bin/bash


NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	BAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/1_SAMtoBAM/MG9-"$NUM".bam"
	BED=/media/cc/A/Alicia/Genome_files/Josemi/ncbiRefSeqCurated.bed
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/2_Intersect_RefSeqBED_BAM/MG9-"$NUM".bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done


