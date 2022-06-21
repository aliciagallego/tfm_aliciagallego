#!/bin/bash

# This script merges the list of genes with significant expression in the 4sU-DRBseq experiment (previously computed) with a reference transcriptome
# and saves the longest version of each transcript including a number of parameters needed for the elongation rate calculation, 
# i.e.: ID, chr, strand, start coord, end coord, CDS start, CDS end, Exon number, Exon starts, Exon ends, splice variants
# The generated file is the Input_gene_list.txt used as input for elongation rate calculation

OUT="/path/Input_genes"

paste /path/Intersect_RefSeq20Kb_Normalized_MG9_12.txt /path/Intersect_RefSeq20Kb_Normalized_MG9_14.txt /path/Intersect_RefSeq20Kb_Normalized_MG9_16.txt /path/Intersect_RefSeq20Kb_Normalized_MG9_18.txt | awk -v OFS="\t" '($6>0.7 && $12>0.7 && $18>0.7 && $24>0.7 && $3-$2>20000){print $1,$2,$3,$4,1,$5}' > $OUT/MYGENES.bed

cat $OUT/MYGENES.bed | awk '{print $4}' | sort | uniq | while read LINE
do
	printf "chrZ\t1\t1\t1\t1\t1\t1\n" > REMOVE.ME
	cat /path/RefSeqGenes_mm10/ncbiRefSeqCurated.txt | grep -P "\t$LINE\t" | while read LINE2
	do
		PREVLENGTH=$(cat REMOVE.ME | head -1 | awk '{print $6-$5}')
		LENGTH=$(echo $LINE2 | awk '{print $6-$5}')
		if [ $LENGTH -gt $PREVLENGTH ]
		then
			echo $LINE2 | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > REMOVE.ME
		fi
	done
	cat REMOVE.ME >> REMOVE.ME2
done

printf "RN\tID\tChromosome\tOrientation\tStart\tEnd\tCDS_start\tCDS_end\tExon_number\tExon_starts\tExon_ends\tZ\tName\tcmpl1\tcmpl2\tSP\n" > $OUT/Input_genes_20Kb.txt
cat REMOVE.ME2 | grep -v chrZ | sort -k3,3 -k5,5n -k6,6n | sed 's/chrX/20/g' | sed 's/chrY/21/g' | sed 's/chr//g' | sed "s/\t+\t/\t1\t/g" | sed "s/\t-\t/\t2\t/g" >> $OUT/Input_genes_20Kb.txt
rm -f REMOVE.ME* MYGENES.bed
