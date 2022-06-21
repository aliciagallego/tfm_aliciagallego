#!/bin/bash


for NUM in TKO_I TKO_II TKO_III WT_I WT_II WT_III
do
	BAM="/media/cc/A/Josemi/NGS/cheRNA/BAMs/"$NUM".bam"
	BED=/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed
	SALIDA="/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/1_Intersect_RefSeq/"$NUM"_RefSeq_intersect.bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done


