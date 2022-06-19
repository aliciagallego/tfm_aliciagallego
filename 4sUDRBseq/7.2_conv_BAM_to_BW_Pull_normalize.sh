#!/bin/bash

# This script generates BigWig files from bam files 
# -bs: size of the bins, in bases, for the output of the bigwig file (default: 50)

#NUMS=$(TKO0 TKO5 WT0 WT5);

for NUM in TKO0 TKO5 WT0 WT5
do
        BAM="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/Pull/"$NUM"_*.bam"

	FILE=$(basename $BAM)
	NAME=${FILE%.*}

	OUTFW="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/Pull/"$NAME"_FW.bw"
	OUTRV="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/Pull/"$NAME"_RV.bw"

	#samtools index $BAM
	#bamCoverage -bs 1 -b $BAM -o $OUTPUT
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand forward -o $OUTFW
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand reverse -o $OUTRV
done
wait

