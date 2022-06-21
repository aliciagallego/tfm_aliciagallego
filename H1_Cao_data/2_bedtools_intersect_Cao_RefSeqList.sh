#!/bin/bash

AFILE=/path/H1_Cao/Intersect/Input_genes_RefSeq_Long_List_2KbTSS.bed
BAMS=/path/H1_Cao/Alignment/Pull/

for FILE in $BAMS/*.bam
do

	NAME=$(basename $FILE | sed 's/_s.*//') 
       	bedtools intersect -a $AFILE -b $FILE -c > path/H1_Cao/Intersect/$NAME'_2kbTSS'.bed

done
