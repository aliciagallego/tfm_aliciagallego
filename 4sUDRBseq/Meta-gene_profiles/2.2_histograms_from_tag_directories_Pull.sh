#!/bin/bash

# Tag directory: is a directory that contains several files describing your data. 
# The input file for making a tag directory is the output of the mapping process *.sam

for NUM in 11-13 12-14 15-17 16-18
do

	TAG="/path/4suDRB/Pull/TagDirectories_Tag/MG9-"$NUM"_TagDirectory/"

	OUT="/path/4suDRB/Pull/Histograms//MG9-"$NUM"_Histogram"      

	annotatePeaks.pl /path/4suDRB/RefSeq_LongList_TSS_2Kb50Kb.bed none -size 52000 -hist 20 -ghist -d $TAG > $OUT/"MG9-"$NUM"_output.txt"     	 	


done
wait
