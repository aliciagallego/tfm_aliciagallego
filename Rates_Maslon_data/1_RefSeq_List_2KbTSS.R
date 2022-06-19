#!/usr/bin/env Rscript

#############################################################
## Include 2Kb region at TSS (start) and TTS (end) extremes #
#############################################################

# -------
# Paths |
# -------
ensembl_path <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_data/Ensembl_mm10_transcripts.bed"
rates_out <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_output/1_Prepare_Ensembl_list/"

# -----------
# Open data |
# -----------
ensembl <- read.table(ensembl_path,h=F,sep="\t",stringsAsFactors=FALSE)

# -------------------------------------
# Transform start - 2Kb and end + 2Kb |
# -------------------------------------
count <- 0
for(i in 1:nrow(ensembl)) {
  if(ensembl[i,6] == '+' | ensembl[i,6] == '-') {
    ensembl[i,2]=(ensembl[i,2]-2000)
    ensembl[i,3]=(ensembl[i,3]+2000)
  }
  else {
    ensembl[i,8]=NA
    ensembl[i,9]=NA
  }
  count = count+1
}
print(count)

# --------------------------------------------------------------------------------------
# Remove transcripts with negative coords (these are mitochondrial and rare sequences) |
# --------------------------------------------------------------------------------------
aa<- ensembl$V1=="chrMT"
bb <- ensembl$V1!="chrMT"
bb <- ensembl[!ensembl$V1 == "chrMT",]
ensembl2 <- ensembl[!ensembl$V2 < 0,]

# -----------
# Save data |
# -----------
write.table(ensembl2, 
            file = paste0(rates_out,"Ensembl_mm10_transcripts_2Kb.bed"),
            quote = F, 
            sep="\t", 
            col.names = F, 
            row.names = F)
