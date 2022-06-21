#!/bin/bash

# This script computes coverage per base for each chromosome
# Inputs: SAM alignments and a list of chromosome sizes
# Outputs: profile_$chr 

process () {

	cat $1 | grep -vP "^@" | awk -v OFS="\t" '($5>5 && $7=="=" && $9>0 && $9<800){print $3,$4,$4+$9}' | sort -k1,1 -k2,2n -k3,3n > REMOVE.ME.bed
	NUMS=$(seq 1 19); 
	for NUM in $NUMS X Y
        do
                cat REMOVE.ME.bed | grep -P "chr$NUM\t" > "chr"$NUM".bed"
                cat /path/Genome_files/mm10/chrom.sizes | grep -P "chr$NUM\t" > GENOME$NUM
		#{printf ("%-20d\n", $3)} removes scientific notation for big numbers
                genomeCoverageBed -d -i "chr"$NUM".bed" -g GENOME$NUM | awk '{printf ("%-20d\n", $3)}' | tr "\n" "\t" > profile_$NUM
                rm -f GENOME$NUM
        done
}

THISDIR=/path/Rate_calculation

for FILE in /path/Pull/*.sam
do
	NAME=$(basename $FILE .sam)
	mkdir $NAME
	cd $NAME
	process $FILE &
	cd $THISDIR
done
wait

for DIR in MG*
do
	mv $DIR/profile_X $DIR/profile_20
	mv $DIR/profile_Y $DIR/profile_21
done
