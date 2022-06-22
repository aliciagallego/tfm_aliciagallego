#!/bin/bash

AFILE=/path/RNAPII/Gene_selection_RefSeq/Input_genes_RefSeq_Long_List_Promoters.bed
BAM=/path/RNAPII/Alignment

for FILE in $BAM/*.bam
do
	NAME=$(basename $FILE | sed 's/_R.*//')
       	bedtools intersect -c -wa -a $AFILE -b $FILE  > /path/RNAPII/Intersect_bedtools/Promoters/$NAME'_wa_PR_500bp'.bed

done
