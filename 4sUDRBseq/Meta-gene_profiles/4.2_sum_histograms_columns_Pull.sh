#!/bin/bash

for NUM in 11-13 12-14 15-17 16-18

do
	FRAC=0.01	
	INPUT="/path/4suDRB/Pull/Histograms/MG9-"$NUM"_Histogram/MG9-"$NUM"_matrix_"$FRAC".txt"
	OUT="/path/4suDRB/Pull/Histograms/MG9-"$NUM"_Histogram"
	OUT2="/path/4suDRB/Pull/Histograms/"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/"MG9-"$NUM"_matrix_sum_"$FRAC".txt"
	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'>> $OUT2/"sum_all_matrix_sum_"$FRAC".txt"

done
wait
