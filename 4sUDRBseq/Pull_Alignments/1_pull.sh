#!/bin/bash


process () {
	gzip -cd $1 > $3
	gzip -cd $2 >> $3
	gzip $3
}


for NUM in 11 12 15 16
do
	FILE1="/media/cc/B/Josemi/TTseq_Feb2022/Reads_project_MG10/MG10_RNA_DRB-4sU/MG10_RNA/*S"$NUM"_S"$NUM"_R1_001.fastq.gz"
	FILE2="/media/cc/B/Josemi/TTseq_Feb2022/Reads_project_MG10/MG10_RNA_DRB-4sU/MG10_RNA/*S"$NUM"_S"$NUM"_R2_001.fastq.gz"
	
	NUMPAIR=$(($NUM + 2))
	PAIR1="/media/cc/B/Josemi/TTseq_Feb2022/Reads_project_MG10/MG10_RNA_DRB-4sU/MG10_RNA/*S"$NUMPAIR"_S"$NUMPAIR"_R1_001.fastq.gz"
        PAIR2="/media/cc/B/Josemi/TTseq_Feb2022/Reads_project_MG10/MG10_RNA_DRB-4sU/MG10_RNA/*S"$NUMPAIR"_S"$NUMPAIR"_R2_001.fastq.gz"

	SALIDA1=$(basename $FILE1 .gz)
	SALIDA2=$(basename $FILE2 .gz)

	#process $FILE1 $PAIR1 $SALIDA1 &
	#process $FILE2 $PAIR2 $SALIDA2 &

	echo $FILE1
	echo $FILE2
done
wait
