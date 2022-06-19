#!/usr/bin/env Rscript

#############################################################
## Include 2Kb region at TSS (start) and TTS (end) extremes #
#############################################################

# -------
# Paths |
# -------
refseq_path <- "/home/agj/Documentos/TFM/RNApolII_2/RefSeq_genes.bed"
output <- "/home/agj/Documentos/TFM/meRIP_2/meRIP2_output/1_Gene_selection_RefSeq/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE)

# -------------------------------------
# Transform start - 2Kb and end + 2Kb |
# -------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+' | refseq[i,6] == '-') {
    refseq[i,2]=(refseq[i,2]-2000)
    refseq[i,3]=(refseq[i,3]+2000)
  }
  else {
    refseq[i,8]=NA
    refseq[i,9]=NA
  }
  count = count+1
}
print(count)

# --------------------------------------------------------------------------------------
# Remove transcripts with negative coords (these are mitochondrial and rare sequences) |
# --------------------------------------------------------------------------------------
aa<- refseq$V1=="chrMT" #0
refseq <- refseq[!refseq$V2 < 0,] #2 genes from rare chr

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"RefSeq_LongList_2Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)