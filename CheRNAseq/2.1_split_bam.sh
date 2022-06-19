#!/bin/bash

for FILE in TKO WT

do
	BAM="/media/cc/A/Josemi/NGS/cheRNA/BAMs/Pulls/"$FILE".bam"
	SALIDA="/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/1_Intersect_RefSeq/Splits/"$NAME
	NAMEE=$(basename $BAM)
	NAME=${NAMEE%.*}"_Pull"
	
	LINES=$(samtools view -c $BAM)  # get number of lines
	HALF=$(($LINES/2)) # get half number of lines

	samtools view -H $BAM > $NAME"_split1.sam" # get headers and create split 1
	samtools view $BAM | head -n $HALF >>  $NAME"_split1.sam"  # paste first half of the lines below the headers
	samtools view -S -b $NAME"_split1.sam" > $NAME"_split1.bam" # from sam to bam
	rm $NAME"_split1.sam" # remove sam


	samtools view -H $BAM > $NAME"_split2.sam" # get headers and create split 2
	samtools view $BAM | tail -n $HALF >>  $NAME"_split2.sam"  # paste second half of the lines below the headers
	samtools view -S -b $NAME"_split2.sam" > $NAME"_split2.bam" # from sam to bam
	rm $NAME"_split2.sam" # remove sam

done

