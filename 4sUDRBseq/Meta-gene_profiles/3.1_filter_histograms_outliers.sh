#!/bin/bash

NUMS=$(seq 11 18);

for NUM in $NUMS
do
	INPUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5_Histograms2/ghist/MG9-"$NUM"_Histogram/MG9-"$NUM"_output.txt"
	OUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/5_Histograms2/ghist/MG9-"$NUM"_Histogram"
	FRAC=0.001	

	# obtener el numero total de filas (genes) sin contar la fila del encabezado
	MYNUM=$(cat $INPUT | sed '1d' | wc -l)

	# multiplicar el numero total de filas por 0.05 y sumar 1
	# bc -l indica usar formulas matemáticas
	# sed quita los decimales 4.55 -> 4
	MYNUM=$(echo "($MYNUM * $FRAC) + 1" | bc -l | sed 's/\..*//')

	# Threshold
	# sed "1d" Prints the contents of file excluding the first line to the standard output
	# sumar todas las columnas por fila
	# sed 's/\..*//':  quita los decimales
	# sort -nr: ordena de mayor a menor
	# head -n $MYNUM: selecciona los $MYNUM primeros
	# tail -1: selecciona el último de los valores $MYNUM (el número a partir del cual se van a eliminar los valores altos o threshold)

	THRES=$(cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print sum; sum=0}' | sed 's/\..*//' | sort -nr  | head -n $MYNUM | tail -1)

	printf "up on the threshold\t"$MYNUM"\n" > $OUT/"MG9-"$NUM"_summary_0.001.txt"
	printf "threshold\t"$THRES"\n" >> $OUT/"MG9-"$NUM"_summary_0.001.txt"

	# imprimir los encabezados
	cat $INPUT | head -1 > $OUT/"MG9-"$NUM"_matrix_0.001.txt"

	# sed "1d": Prints the contents of file excluding the first line to the standard output
	# sumar todas las columnas por fila
	# seleccionar aquellas filas cuya columna sumatoria es menor al threhold anteriormente calculado
	cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print $0,sum; sum=0}' | awk -v OFS="\t" -v thres=$THRES '($NF<thres){printf $1"\t"; for (i=2; i<NF; i++) printf $i"\t"; printf "\n"}' >> $OUT/"MG9-"$NUM"_matrix_0.001.txt"     	 	

done
wait
