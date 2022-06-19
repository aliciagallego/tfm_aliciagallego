# Tag directory: is a directory that contains several files describing your data. In this case the input file for making a tag directory are the .sorted.bam files (obtained from *.sam)

# makeTagDirectory is a function from HOMER

for NUM in 11-13 12-14 15-17 16-18
do	
	INPUT="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/2_Sorting/Pull/*_MG9-"$NUM".bam"
	
	OUTPUTDIR="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/4.2_TagDirectories/MG9-"$NUM"_TagDirectory/"        	

	# dudo si debo usar los parametros -flip -sspe 
	makeTagDirectory $OUTPUTDIR -flip -sspe $INPUT -format sam 
done
wait


