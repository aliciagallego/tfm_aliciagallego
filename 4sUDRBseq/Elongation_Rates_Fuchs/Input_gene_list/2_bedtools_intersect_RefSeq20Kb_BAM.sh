#!/bin/bash

NUMS=$(seq 11 18); 

for NUM in 12 14 16 18
do
	BAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/6_SAM_filtered_q1/MG9-"$NUM"_filtered.bam"
	BED=/media/cc/B/Josemi/TTseq_Feb2022/TTseq_data/RefSeq_LongList_TSS_20Kb.bed
	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/2_Intersect_RefSeqBED20Kb_BAM/MG9-"$NUM".bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done
