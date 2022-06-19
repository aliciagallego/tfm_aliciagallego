#!/bin/bash


for NUM in TKO_Pull WT_Pull
do
	BAM1="/media/cc/A/Josemi/NGS/cheRNA/BAMs/Pulls/"$NUM"_split1.bam"
	BAM2="/media/cc/A/Josemi/NGS/cheRNA/BAMs/Pulls/"$NUM"_split2.bam"
	BED=/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed
	SALIDA1="/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/1_Intersect_RefSeq/"$NUM"_RefSeq_intersect1.bed"
	SALIDA2="/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/1_Intersect_RefSeq/"$NUM"_RefSeq_intersect2.bed"

	bedtools intersect -a $BED -b $BAM1 -c -s > $SALIDA1
	bedtools intersect -a $BED -b $BAM2 -c -s > $SALIDA2
done


