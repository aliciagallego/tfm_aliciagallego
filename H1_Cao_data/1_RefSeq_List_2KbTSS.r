#!/usr/bin/env Rscript

####################################################################################
## Get promoters coords (+-2 Kp) from RefSeq transcript start/end coords Long List #
####################################################################################

# -------
# Paths |
# -------
refseq_path <- "/path/Genome_files/RefSeq_genes.bed"
output <- "/path/H1_Cao/Intersect/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=F)
nrow(refseq)

# --------------------------------------------------
# Generate new start/end coords from previous ones |
# --------------------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+') {
    refseq[i,7]=(refseq[i,2]-2000)
    refseq[i,8]=(refseq[i,2]+2000)
  }
  else if(refseq[i,6] == '-'){
    refseq[i,7]=(refseq[i,3]-2000)
    refseq[i,8]=(refseq[i,3]+2000)
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

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"Input_genes_RefSeq_Long_List_2KbTSS.bed"),
            quote = F, sep="\t", col.names = F, row.names = F) 
