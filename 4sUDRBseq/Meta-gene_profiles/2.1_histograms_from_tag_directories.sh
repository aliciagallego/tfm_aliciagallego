# Tag directory: is a directory that contains several files describing your data. 
# The input file for making a tag directory is the output of the mapping process *.sam


NUMS=$(seq 11 18);

for NUM in $NUMS
do
	TAG="/path/4suDRB/TagDirectories/MG9-"$NUM"_TagDirectory/"
	OUT="/path/4suDRB/Histograms/MG9-"$NUM"_Histogram"      

	annotatePeaks.pl /path/4suDRB/RefSeq_LongList_TSS_2Kb50Kb.bed none -size 52000 -hist 20 -ghist -d $TAG > $OUT/$3"_output.txt"     	 	


done
wait
