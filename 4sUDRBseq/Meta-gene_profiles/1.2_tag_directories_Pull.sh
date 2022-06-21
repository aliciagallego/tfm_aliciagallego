# Tag directory: is a directory that contains several files describing your data. 
# In this case the input file for making a tag directory are the .sorted.bam files (obtained from *.sam)
# makeTagDirectory is a function from HOMER

for NUM in 11-13 12-14 15-17 16-18
do	
	INPUT="/path/4suDRB/Pull/*_MG9-"$NUM".bam"
	
	OUTPUTDIR="/path/4suDRB/Pull/TagDirectories/MG9-"$NUM"_TagDirectory/"        	
	
	makeTagDirectory $OUTPUTDIR -flip -sspe $INPUT -format sam 
done
wait
