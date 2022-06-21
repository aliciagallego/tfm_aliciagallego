#!/bin/bash

# This script generates BigWig files from bam files 
# -bs: size of the bins, in bases, for the output of the bigwig file (default: 50)
# --normalizeUsing RPKM normalization is used

NUMS=$(seq 11 18);

for NUM in $NUMS
do
        BAM="/path/4sUDRB/MG9-"$NUM"*.bam"

	FILE=$(basename $BAM)
	NAME=${FILE%.*}

	OUTFW="/path/4sUDRB/"$NAME"_Watson.bw"
	OUTRV="/path/4sUDRB/"$NAME"_Crick.bw"
	
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand forward -o $OUTFW
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand reverse -o $OUTRV
done
wait
