#!/usr/bin/env Rscript

###############################################################################
## Get gene body coords (+500bp) from RefSeq transcript start coord Long List #
###############################################################################

# -------
# Paths |
# -------
refseq_path <- "/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed"
output <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Gene_selection_RefSeq/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=F)

# --------------------------------------------------
# Generate new start/end coords from previous ones |
# --------------------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,6] == '+') {
    refseq[i,2]=(refseq[i,2]+500)
  }
  else if(refseq[i,6] == '-'){
    refseq[i,3]=(refseq[i,3]-500)
  }
  count = count+1
}
print(count)

# -------------------------------------
# Remove genes less than 500 bp  long |
# -------------------------------------
refseq <- subset(refseq, refseq[,3]-refseq[,2] > 0)

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"Input_genes_RefSeq_Long_List_GeneBody.bed"),
            quote = F, sep="\t", col.names = F, row.names = F) 
