#!/usr/bin/env Rscript

###################################################
## Normalize 4sUDRB expression (RefSeq gene list) # 
###################################################

# The objective of this script is the normalization by total number of reads per experiment of the transcriptome RefSeq list based on 4sUDRBseq reads
# It also adds the correct star and end coords for each gene (since they were transformed to just consider the first 20 Kb)

# -----------
# Libraries |
# -----------
library(dplyr)

# -------
# Paths |
# -------
refseq_path <- "/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed"
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/2_Intersect_RefSeqBED20Kb_BAM/"
output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/3_Normalization/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("Chr","Start","End","Gene_name","NA1","Strand"))

TTseq_list = list.files(TTseq_path, pattern="*.bed")

for (i in seq_along(TTseq_list)) {
  filename <- sub(".bed", "", TTseq_list[i])
  filename <- sub("-", "_", filename)
  df <- read.table(paste0(TTseq_path,TTseq_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",paste0("TTseq_",filename)))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(filename, df)
}

# ------------------------------------
# Total reads (data from experiment) | Only 5 min samples are taken into account
# ------------------------------------
MG9_12_reads <- 150557906
MG9_14_reads <- 131822641
MG9_16_reads <- 124134871
MG9_18_reads <- 85923486

# ------------------
# 1. Normalization | only by total reads (no transcript size)
# ------------------
MG9_12$TTseq_MG9_12 = (MG9_12$TTseq_MG9_12/MG9_12_reads) * 10000000
MG9_14$TTseq_MG9_14 = (MG9_14$TTseq_MG9_14/MG9_14_reads) * 10000000
MG9_16$TTseq_MG9_16 = (MG9_16$TTseq_MG9_16/MG9_16_reads) * 10000000
MG9_18$TTseq_MG9_18 = (MG9_18$TTseq_MG9_18/MG9_18_reads) * 10000000

# --------------------------------------------------------------------------------------
# Remove transcripts with negative coords (these are mitochondrial and rare sequences) |
# --------------------------------------------------------------------------------------
aa<- refseq$V1=="chrMT" #0
refseq <- refseq[!(refseq$Strand == '-' & (refseq$End-20000<0)),]

# --------------------------------------------
# Recover start and end original coordinates |
# --------------------------------------------
MG9_12_merged<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                       list(MG9_12[,c("Gene_name","Chr","Strand","TTseq_MG9_12")],
                            refseq[,c("Gene_name","Start", "End")]))
MG9_14_merged<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                       list(MG9_14[,c("Gene_name","Chr","Strand","TTseq_MG9_14")],
                            refseq[,c("Gene_name","Start", "End")]))
MG9_16_merged<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                       list(MG9_16[,c("Gene_name","Chr","Strand","TTseq_MG9_16")],
                            refseq[,c("Gene_name","Start", "End")]))
MG9_18_merged<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                       list(MG9_18[,c("Gene_name","Chr","Strand","TTseq_MG9_18")],
                            refseq[,c("Gene_name","Start", "End")]))

MG9_12_merged <- select(MG9_12_merged, Chr, Start, End, Gene_name,Strand,TTseq_MG9_12)
MG9_14_merged <- select(MG9_14_merged, Chr, Start, End, Gene_name,Strand,TTseq_MG9_14)
MG9_16_merged <- select(MG9_16_merged, Chr, Start, End, Gene_name,Strand,TTseq_MG9_16)
MG9_18_merged <- select(MG9_18_merged, Chr, Start, End, Gene_name,Strand,TTseq_MG9_18)

# -----------
# Save data |
# -----------
write.table(MG9_12_merged, 
            file = paste0(output,"With_real_transcript_size/Intersect_RefSeq20Kb_Normalized_MG9_12_realsize.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_14_merged, 
            file = paste0(output,"With_real_transcript_size/Intersect_RefSeq20Kb_Normalized_MG9_14_realsize.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_16_merged, 
            file = paste0(output,"With_real_transcript_size/Intersect_RefSeq20Kb_Normalized_MG9_16_realsize.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_18_merged, 
            file = paste0(output,"With_real_transcript_size/Intersect_RefSeq20Kb_Normalized_MG9_18_realsize.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
