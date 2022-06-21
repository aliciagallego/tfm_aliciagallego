#!/bin/bash

# bowtie2 alignment of fastq reads and mm10 indexed genome
# -x basename of the index for the reference genome
# -1 pair 1
# -2 pair 2
# -N number of mismatches to allowed in a seed alignment during multiseed alignment (0 or 1)

NUMS=$(seq 11 18); 

for NUM in $NUMS
do
	READS1="/path/4suDRB/*S"$NUM"_R1_001.fastq.gz"
	READS2="/path/4suDRB/*S"$NUM"_R2_001.fastq.gz"

	SALIDA="/path/4suDRB/Alignments/MG9-"$NUM".sam"

	bowtie2 -x /path/Genome_files/mm10/Bowtie_index -1 $READS1 -2 $READS2 -S $SALIDA -N 1
done
wait
