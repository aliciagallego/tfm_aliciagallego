#!/bin/bash

paste /media/cc/A/Josemi/PRUEBA_metlab/Transcriptome_FPKMs/MG9-3_RPKMs.bed /media/cc/A/Josemi/PRUEBA_metlab/Transcriptome_FPKMs/MG9-8_RPKMs.bed /media/cc/A/Josemi/PRUEBA_metlab/Transcriptome_FPKMs/MG9-13_RPKMs.bed /media/cc/A/Josemi/PRUEBA_metlab/Transcriptome_FPKMs/MG9-18_RPKMs.bed | grep -iv inter | grep -iv intra | grep -iv prev | grep -iv bidi | awk -v OFS="\t" '($5>0.1 && $11>0.1 && $17>0.1 && $23>0.1 && $3-$2>20000){print $1,$2,$3,$4,1,$6}' > MYGENES.bed

cat MYGENES.bed | awk '{print $4}' | sed 's/:.*//' | sort | uniq | while read LINE
do
	printf "chrZ\t1\t1\t1\t1\t1\t1\n" > REMOVE.ME

	cat /home/cc/JoseMiguel/Genome_files/RefSeqGenes_mm10/ncbiRefSeqCurated.txt | grep -P "\t$LINE\t" | while read LINE2
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

printf "RN\tID\tChromosome\tOrientation\tStart\tEnd\tCDS_start\tCDS_end\tExon_number\tExon_starts\tExon_ends\tZ\tName\tcmpl1\tcmpl2\tSP\n" > Input_genes.txt
cat REMOVE.ME2 | grep -v chrZ | sort -k3,3 -k5,5n -k6,6n | sed 's/chrX/20/g' | sed 's/chrY/21/g' | sed 's/chr//g' | sed "s/\t+\t/\t1\t/g" | sed "s/\t-\t/\t2\t/g" >> Input_genes.txt
rm -f REMOVE.ME* MYGENES.bed
