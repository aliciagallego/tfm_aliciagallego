#!/usr/bin/env Rscript

#####################################################################################
## Get promoters coords (+-500bp) from RefSeq transcript start/end coords Long List # 20220214
#####################################################################################

# -------
# Paths |
# -------
refseq_path <- "/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed" # from Josemi (20621 genes)
genebody_path <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Gene_selection_RefSeq/Input_genes_RefSeq_Long_List_GeneBody.bed" # (20546 genes)
output <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Gene_selection_RefSeq/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=F)
genebody <- read.table(genebody_path,h=F,sep="\t",stringsAsFactors=F)
nrow(refseq)
nrow(genebody)
head(genebody)

# --------------------------------------------------
# Generate new start/end coords from previous ones |
# --------------------------------------------------
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
nrow(refseq)

# -----------
# Save data |
# -----------
write.table(refseq, 
            file = paste0(output,"Input_genes_RefSeq_Long_List_Promoters.bed"),
            quote = F, sep="\t", col.names = F, row.names = F) 