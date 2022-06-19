#!/bin/bash


## Transform transcript data in txt to bed and maintain bed format

# tail -n +2: prints all the file starting from file 2


# Three first columns are mandatory in BED files and must be 1(chr), 2(Start), 3(End). Column 4 has free format so we take advantage of this and introduce all data inside. Column 6 must be the strand information (+/-)

# info variable: allows to put all TXT columns accumulated in the fourth coumn of the new BED file and separated by ;; just in case we may need this information later on

# gene name is included in column 7. This is not a BED file but we will need gene name later on

INPUT=/media/cc/A/Alicia/Genome_files/Josemi/ncbiRefSeqCurated.txt
OUTPUT=/media/cc/A/Alicia/Genome_files/Josemi

cat $INPUT | awk '{info=$1; for(i=2;i<=NF;i++){info=info";;"$i}; print $3"\t"$5"\t"$6"\t"info"\tNA\t"$4"\t"$13}' > $OUTPUT/ncbiRefSeqCurated.bed
