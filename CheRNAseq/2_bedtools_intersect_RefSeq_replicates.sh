#!/bin/bash


for NUM in TKO_I TKO_II TKO_III WT_I WT_II WT_III
do
	BAM="/path/cheRNA/Alignments/"$NUM".bam"
	BED=/path/Genome_files/RefSeq_genes.bed
	SALIDA="/path/cheRNA/Intersect_RefSeq/"$NUM"_RefSeq_intersect.bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done
