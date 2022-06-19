#!/usr/bin/env Rscript

############################################
## Retain the first 20 Kb from TSS (start) # 20220325
############################################

# -------
# Paths |
# -------
refseq_path <- "/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed"
output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_data/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("Chr","Start","End","Gene_name","NA1","Strand"))
head(refseq)
# ----------------------------------------
# Transform start - 2Kb and start + 50Kb |
# ----------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+') {
    refseq[i,7]=(refseq[i,2])
    refseq[i,8]=(refseq[i,2]+20000)
  }
  else if(refseq[i,6] == '-'){
    refseq[i,7]=(refseq[i,3]-20000)
    refseq[i,8]=(refseq[i,3])
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
head(refseq)

# --------------------------------------------------------------------------------------
# Remove transcripts with negative coords (these are mitochondrial and rare sequences) |
# --------------------------------------------------------------------------------------
aa<- refseq$V1=="chrMT" #0
aa
refseq <- refseq[!refseq$Start < 0,] #2 genes
nrow(refseq)
# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"RefSeq_LongList_TSS_20Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
