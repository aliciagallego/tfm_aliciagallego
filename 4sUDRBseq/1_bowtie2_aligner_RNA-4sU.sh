#!/bin/bash

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	READS1="/media/cc/B/Josemi/TTseq_Feb2022/Reads_project_MG10/MG10_RNA_DRB-4sU/MG10_RNA/*S"$NUM"_R1_001.fastq.gz"
	READS2="/media/cc/B/Josemi/TTseq_Feb2022/Reads_project_MG10/MG10_RNA_DRB-4sU/MG10_RNA/*S"$NUM"_R2_001.fastq.gz"

	SALIDA="/media/cc/B/Josemi/TTseq_Feb2022/Alignments/MG10_RNA_DRB-4sU/MG9-"$NUM".sam"

	bowtie2 -x /home/cc/JoseMiguel/Genome_files/mm10/Bowtie_index -1 $READS1 -2 $READS2 -S $SALIDA -N 1
done
wait
