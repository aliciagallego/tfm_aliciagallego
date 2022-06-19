#!/bin/bash

AFILE=/media/cc/B/Alicia/H1_Cao/H1_Cao_output/1_bedtools_intersect/Input_genes_RefSeq_Long_List_2KbTSS.bed
BAMS=/media/cc/A/Josemi/Cao_H1/Alignment/Pull/

for FILE in $BAMS/*.bam
do

	NAME=$(basename $FILE | sed 's/_s.*//') 
       	bedtools intersect -a $AFILE -b $FILE -c > /media/cc/B/Alicia/H1_Cao/H1_Cao_output/1_bedtools_intersect/$NAME'_2kbTSS'.bed

done

