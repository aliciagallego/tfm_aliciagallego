#!/bin/bash

AFILE=/media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/1_Gene_selection_RefSeq/RefSeq_LongList_2Kb.bed

INPUTS=/media/cc/A/Alicia/NGS/meRIP/meRIP_data/Alignments/Inputs/*001

for DIR in $INPUTS
do
	NAME=$(basename $DIR | sed 's/_0.*/_2Kb/')
	FILE=$DIR/accepted_hits.bam	
	#OUTDIR=$(basename $DIR | sed 's/_0.*/bedtools/')
	#mkdir $OUTDIR
	
	#bedtools intersect -a $AFILE -b $FILE -c -s > $OUTDIR/$NAME.bed
	bedtools intersect -a $AFILE -b $FILE -c -s > /media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/2_Intersect_RefSeq_BAM/$NAME.bed
	
done

IPS=/media/cc/A/Alicia/NGS/meRIP/meRIP_data/Alignments/IPs/*001

for DIR in $IPS
do
	NAME=$(basename $DIR | sed 's/_0.*/_2Kb/')
	FILE=$DIR/accepted_hits.bam	

	bedtools intersect -a $AFILE -b $FILE -c -s > /media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/2_Intersect_RefSeq_BAM/$NAME.bed
	
done
