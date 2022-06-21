# Tag directory: is a directory that contains several files describing your data. 
# The input file for making a tag directory is the output of the mapping process *.sam
# makeTagDirectory is a function from HOMER

NUMS=$(seq 11 18);

for NUM in $NUMS
do
	INPUT="/path/4suDRB/MG9-"$NUM".sam"
	OUTPUTDIR="/path/4suDRB/TagDirectories/MG9-"$NUM"_TagDirectory/"        	

	makeTagDirectory $OUTPUTDIR -flip -sspe $INPUT -format sam 
done
wait
