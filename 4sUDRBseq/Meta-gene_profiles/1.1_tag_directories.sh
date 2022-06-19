# Tag directory: is a directory that contains several files describing your data. The input file for making a tag directory is the output of the mapping process *.sam

# makeTagDirectory is a function from HOMER

NUMS=$(seq 11 18);

for NUM in $NUMS
do
	INPUT="/media/cc/B/Josemi/TTseq_Feb2022/Alignments/MG10_RNA_DRB-4sU/MG9-"$NUM".sam"
	OUTPUTDIR="/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Visualize_alignments/MG10_RNA_DRB-4sU/4_TagDirectories/MG9-"$NUM"_TagDirectory/"        	

	# dudo si debo usar los parametros -flip -sspe 
	makeTagDirectory $OUTPUTDIR -flip -sspe $INPUT -format sam 
done
wait


