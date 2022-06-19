#!/bin/bash

AFILE=/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Gene_selection_RefSeq/Input_genes_RefSeq_Long_List_GeneBody.bed

NEW=/media/cc/A/Alicia/NGS/RNApolII/RNApolII_data/Archivos_cluster/New_seq
OLD=/media/cc/A/Alicia/NGS/RNApolII/RNApolII_data/Archivos_cluster/Old_seq

for FILE in $NEW/*.bam
do

	NAME=$(basename $FILE | sed 's/_S.*//') 
       	bedtools intersect -c -wa -a $AFILE -b $FILE > /media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Intersect_RefSeq_BAM/1_Intersect_bedtools/GeneBody/$NAME'_wa_GB_500bp'.bed

done

for FILE in $OLD/*.bam
do
	NAME=$(basename $FILE | sed 's/_R.*//')
       	bedtools intersect -c -wa -a $AFILE -b $FILE  > /media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Intersect_RefSeq_BAM/1_Intersect_bedtools/GeneBody/$NAME'_wa_GB_500bp'.bed

done


