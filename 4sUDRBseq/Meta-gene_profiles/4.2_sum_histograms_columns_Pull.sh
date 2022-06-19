#!/bin/bash

for NUM in 11-13 12-14 15-17 16-18

do
	FRAC=0.01	
	INPUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5.2_Histograms_pull/ghist/MG9-"$NUM"_Histogram/MG9-"$NUM"_matrix_"$FRAC".txt"
	OUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5.2_Histograms_pull/ghist/MG9-"$NUM"_Histogram"
	OUT2="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5.2_Histograms_pull/ghist"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/"MG9-"$NUM"_matrix_sum_"$FRAC".txt"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'>> $OUT2/"sum_all_matrix_sum_"$FRAC".txt"

	# imprimir los encabezados
	#cat $INPUT | head -1 > $OUT/"MG9-"$NUM"_matrix_sum.txt"
	
	# tail -n +2 es como sed '1d': coge toda la tabla menos el encabezado
	# awk va recorriendo cada columna y metiendo en un array llamado 'a' la suma. al final de todo imprime los valores de 'a'
	# en vez de printarlos uno a uno los metes todos en una variable y los printa al final
	#tail -n +2 $INPUT | awk '{for(i=2;i<=NF;i++){a[i]+=$i}}END{str=a[2]; for(i=3;i<=NF;i++){str=str"\t"a[i]} print str}' >> $OUT/"MG9-"$NUM"_matrix_sum.txt"  

done
wait

