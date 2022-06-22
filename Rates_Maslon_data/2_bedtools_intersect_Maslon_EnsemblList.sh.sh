#!/bin/bash

AFILE=/path/Maslon2019/Genome_data/Ensembl_mm10_transcripts_2Kb.bed

INPUTS=/path/meRIP/Alignments/Inputs/*001

for DIR in $INPUTS
do
	NAME=$(basename $DIR | sed 's/_0.*/_2Kb/')
	FILE=$DIR/accepted_hits.bam	

	bedtools intersect -a $AFILE -b $FILE -c -s > /path/Maslon2019/Intersect_bedtools/$NAME.bed
done


IPS=/path/meRIP/Alignments/IPs/*001

for DIR in $IPS
do
	NAME=$(basename $DIR | sed 's/_0.*/_2Kb/')
	FILE=$DIR/accepted_hits.bam	

	bedtools intersect -a $AFILE -b $FILE -c -s > /path/Maslon2019/Intersect_bedtools/$NAME.bed
done
