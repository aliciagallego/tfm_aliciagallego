#!/bin/bash


NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	BAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/6_SAM_filtered_-q1/MG9-"$NUM"_filtered.bam"
	BED=/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/2.2_Intersect_RefSeqBEDcurated_BAM/MG9-"$NUM".bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done


