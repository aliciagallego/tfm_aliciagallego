#!/bin/bash

AFILE=/path/RNAPII/Gene_selection_RefSeq/Input_genes_RefSeq_Long_List_GeneBody.bed
BAM=/path/RNAPII/Alignment

for FILE in $BAM/*.bam
do

	NAME=$(basename $FILE | sed 's/_S.*//') 
       	bedtools intersect -c -wa -a $AFILE -b $FILE > /path/RNAPII/Intersect_bedtools/GeneBody/$NAME'_wa_GB_500bp'.bed

done
