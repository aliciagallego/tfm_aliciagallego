#!/usr/bin/env Rscript

#####################################################################################
## Get promoters coords (+-500bp) from RefSeq transcript start/end coords Long List #
#####################################################################################

# -------
# Paths |
# -------
refseq_path <- "/path/Genome_files/RefSeq_genes.bed"
output <- "/path/RNAPII/Gene_selection_RefSeq/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=F)

# ---------------------------------------------------
# Transform start = start-500b and end = start+500b |
# ---------------------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+') {
    refseq[i,7]=(refseq[i,2]-500)
    refseq[i,8]=(refseq[i,2]+500)
  }
  else if(refseq[i,6] == '-'){
    refseq[i,7]=(refseq[i,3]-500)
    refseq[i,8]=(refseq[i,3]+500)
  }
  else {
    refseq[i,7]=NA
    refseq[i,8]=NA
  }
  count = count+1
}
print(count)

# Replace col 2 by col 8 and col 3 by col 9
count <- 0
for(i in 1:nrow(refseq)) {
  refseq[i,2]=refseq[i,7]
  refseq[i,3]=refseq[i,8]
  count = count+1
}
print(count)

# Drop columns 8 and 9
refseq <- refseq[-c(7,8)]

# ---------------------------------------------
# Remove transcripts absent in Gene body list |
# ---------------------------------------------
refseq <- refseq[refseq$V4 %in% genebody$V4, ]

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"Input_genes_RefSeq_Long_List_Promoters.bed"),
            quote = F, sep="\t", col.names = F, row.names = F) 
