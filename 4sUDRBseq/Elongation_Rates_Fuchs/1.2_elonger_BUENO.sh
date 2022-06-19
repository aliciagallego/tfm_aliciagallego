#!/bin/bash

process () {

	cat $1 | grep -vP "^@" | awk -v OFS="\t" '($5>5 && $7=="=" && $9>0 && $9<800){print $3,$4,$4+$9}' | sort -k1,1 -k2,2n -k3,3n > REMOVE.ME.bed
	
	NUMS=$(seq 1 19); 

	for NUM in $NUMS X Y
	#for NUM in {1..19} X Y -> not bash code, maybe python
        do
                cat REMOVE.ME.bed | grep -P "chr$NUM\t" > "chr"$NUM".bed"
                cat /home/cc/JoseMiguel/Genome_files/mm10/chrom.sizes | grep -P "chr$NUM\t" > GENOME$NUM
		#{printf ("%-8d\n", $3)} removes scientific notation for big numbers
                genomeCoverageBed -d -i "chr"$NUM".bed" -g GENOME$NUM | awk '{printf ("%-20d\n", $3)}' | tr "\n" "\t" > profile_$NUM
                rm -f GENOME$NUM
		#rm -f REMOVE.ME.bed
        done
}


THISDIR=/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/1.2_Rate_calculation

for FILE in /media/cc/B/Josemi/TTseq_Feb2022/Alignments/MG10_RNA_DRB-4sU/*sam
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
