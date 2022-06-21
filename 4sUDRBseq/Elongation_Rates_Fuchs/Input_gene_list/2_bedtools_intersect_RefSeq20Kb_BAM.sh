#!/bin/bash

for NUM in 12 14 16 18
do
	BAM="/path/4sUDRB/MG9-"$NUM"_filtered.bam"
	BED=/path/RefSeq_LongList_TSS_20Kb.bed
	SALIDA="/path/Intersect_RefSeqBED20Kb_BAM/MG9-"$NUM".bed"

	bedtools intersect -a $BED -b $BAM -c -s > $SALIDA
done
