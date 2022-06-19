#!/bin/bash

for NUM in 11-13 12-14 15-17 16-18
do

	TAG="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/4.2_TagDirectories/MG9-"$NUM"_TagDirectory/"

	OUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5.2_Histograms_pull/ghist/MG9-"$NUM"_Histogram"      

	annotatePeaks.pl /media/cc/B/Josemi/TTseq_Feb2022/TTseq_data/RefSeq_LongList_TSS_2Kb50Kb_OK.bed none -size 52000 -hist 20 -ghist -d $TAG > $OUT/"MG9-"$NUM"_output.txt"     	 	


done
wait
