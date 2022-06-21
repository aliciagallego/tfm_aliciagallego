#!/usr/bin/env Rscript

###################################################
## Include -2Kb and +50Kb region from TSS (start) #
###################################################

# -------
# Paths |
# -------
refseq_path <- "/path/Genome_files/RefSeq_genes.bed"
output <- "/path/4suDRB/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("Chr","Start","End","Gene_name","NA1","Strand"))

# ----------------------------------------
# Transform start - 2Kb and start + 50Kb |
# ----------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+') {
    refseq[i,7]=(refseq[i,2]-2000)
    refseq[i,8]=(refseq[i,2]+50000)
  }
  else if(refseq[i,6] == '-'){
    refseq[i,7]=(refseq[i,3]-50000)
    refseq[i,8]=(refseq[i,3]+2000)
  }
  else {
    refseq[i,7]=NA
    refseq[i,8]=NA
  }
  count = count+1
}
print(count)

# Replace col 2 by col 7 and col 3 by col 8
count <- 0
for(i in 1:nrow(refseq)) {
  refseq[i,2]=refseq[i,7]
  refseq[i,3]=refseq[i,8]
  count = count+1
}
print(count)

# Drop columns 7 and 8
refseq <- refseq[-c(7,8)]

# --------------------------------------------------------------------------------------
# Remove transcripts with negative coords (these are mitochondrial and rare sequences) |
# --------------------------------------------------------------------------------------
aa<- refseq$V1=="chrMT" #0
refseq <- refseq[!refseq$Start < 0,] #5 genes

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"RefSeq_LongList_TSS_2Kb50Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
