#!/bin/bash

# This script generates BigWig files from bam files 
# -bs: size of the bins, in bases, for the output of the bigwig file (default: 50)
# --normalizeUsing RPKM  normalization is used

for NUM in TKO0 TKO5 WT0 WT5
do
        BAM="/path/4suDRB/"$NUM"_*.bam"

	FILE=$(basename $BAM)
	NAME=${FILE%.*}

	OUTFW="/path/4suDRB/"$NAME"_FW.bw"
	OUTRV="/path/4suDRB/"$NAME"_RV.bw"

	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand forward -o $OUTFW
	bamCoverage -bs 1 -b $BAM --normalizeUsing RPKM --filterRNAstrand reverse -o $OUTRV
done
wait
