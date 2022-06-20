#!/bin/bash

# This script generates BigWig files from bam files 
# -bs: size of the bins, in bases, for the output of the bigwig file (default: 50)

NUMS=$(seq 11 18);

for NUM in $NUMS
do
        BAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/MG9-"$NUM"*.bam"

	FILE=$(basename $BAM)
	NAME=${FILE%.*}

	OUTFW="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/"$NAME"_Watson.bw"
	OUTRV="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/"$NAME"_Crick.bw"
	
	#samtools index $BAM
	#bamCoverage -bs 1 -b $BAM -o $OUTPUT
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand forward -o $OUTFW
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand reverse -o $OUTRV
done
wait
