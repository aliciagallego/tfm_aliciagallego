#!/bin/bash

NUMS=$(seq 11 18);

for NUM in $NUMS
do
	INPUT="/path/4suDRB/Histograms/MG9-"$NUM"_Histogram/MG9-"$NUM"_matrix_0.001.txt"
	OUT="/path/4suDRB/Histograms/MG9-"$NUM"_Histogram"
	OUT2="/path/4suDRB/Histograms/"

	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'> $OUT/"MG9-"$NUM"_matrix_sum_0.001.txt"
	cat $INPUT | awk '{if(NR==1){for(i=2;i<=NF;i++){h[i]=$i} next;} for(i=2;i<=NF;i++){a[i]+=$i}}END{for(i=2;i<=NF;i++){print h[i]"\t"a[i]}}'>> $OUT2/"sum_all_matrix_sum_0.001.txt"

done
wait
