#!/bin/bash

NUMS=$(seq 11 18);

for NUM in $NUMS
do
	INPUT="/path/4suDRB/Histograms/MG9-"$NUM"_Histogram/MG9-"$NUM"_output.txt"
	OUT="/path/4suDRB/Histograms/MG9-"$NUM"_Histogram"
	
	# FRAC value is defined after prior analyses 
	# In this case the 0.1% of the most expressed genes will we removed
	FRAC=0.001	

	# get total number of rows (genes) excluding the header 
	MYNUM=$(cat $INPUT | sed '1d' | wc -l)

	# multiply total row number by 0.05 and sum +1
	# bc -l indicates use math operations
	# sed remoce decimals e.g. 4.55 -> 4
	MYNUM=$(echo "($MYNUM * $FRAC) + 1" | bc -l | sed 's/\..*//')

	# Threshold
	# sed "1d" Prints the contents of file excluding the first line to the standard output
	# sum all columns by row
	# sed 's/\..*//':  removes decimals
	# sort -nr: sorts from higher to lower values
	# head -n $MYNUM: selects the first $MYNUM values
	# tail -1: selects the last value from the $MYNUM selectin (it is the value from which the higher values will be removed or threshold)

	THRES=$(cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print sum; sum=0}' | sed 's/\..*//' | sort -nr  | head -n $MYNUM | tail -1)

	printf "up on the threshold\t"$MYNUM"\n" > $OUT/"MG9-"$NUM"_summary_0.001.txt"
	printf "threshold\t"$THRES"\n" >> $OUT/"MG9-"$NUM"_summary_0.001.txt"

	# prints headers
	cat $INPUT | head -1 > $OUT/"MG9-"$NUM"_matrix_0.001.txt"

	# sed "1d": Prints the contents of file excluding the first line to the standard output
	# sum all columns by row
	# selects all rows which summatory column is lower than the previously computed threshold
	cat $INPUT | sed '1d' | awk -v OFS="\t" '{for(i=2;i<=NF;i++) sum+=$i; print $0,sum; sum=0}' | awk -v OFS="\t" -v thres=$THRES '($NF<thres){printf $1"\t"; for (i=2; i<NF; i++) printf $i"\t"; printf "\n"}' >> $OUT/"MG9-"$NUM"_matrix_0.001.txt"     	 	

done
wait
