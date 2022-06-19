# Tag directory: is a directory that contains several files describing your data. The input file for making a tag directory is the output of the mapping process *.sam


NUMS=$(seq 11 18);

for NUM in $NUMS
do
	TAG="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/4_TagDirectories/MG9-"$NUM"_TagDirectory/"
	OUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5_Histograms2/ghist/MG9-"$NUM"_Histogram"      
	
	# cuidado, ha cambiado el path de la RefSeq_LongList_TSS
	annotatePeaks.pl /media/cc/A/Alicia/NGS/TTseq_Feb2022/TTseq_data/RefSeq_LongList_TSS_2Kb50Kb_OK.bed none -size 52000 -hist 20 -ghist -d $TAG > $OUT/$3"_output.txt"     	 	


done
wait
