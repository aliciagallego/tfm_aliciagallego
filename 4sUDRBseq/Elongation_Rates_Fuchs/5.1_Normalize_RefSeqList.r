#!/usr/bin/env Rscript

#############################
## Normalize (RefSeq 34891) # 20220311
#############################

# The objective of this script is the curation of the transcriptome RefSeq list based on TTseq_Feb2022 alignments
# First: normalization by total number of reads per experiment and transcript size
# Second: filtering by transcript size (discard small transcripts) and expression levels (discard 0 -or low- expressed genes in our data)

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/2.2_Intersect_RefSeqBEDcurated_BAM/"
output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/3_Normalization/"

refseq_path <- "/media/cc/A/Alicia/Genome_files/Josemi/ncbiRefSeqCurated.txt"
refseqbed_path <- "/media/cc/A/Alicia/Genome_files/Josemi/ncbiRefSeqCurated.bed"
refseqshort_path <- "/media/cc/A/Alicia/Genome_files/Josemi/RefSeq_genes.bed"
output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/3_Normalization/"

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("NA1","ID","Chr","Strand","Start","End", "CDS_Start", "CDS_End",
                                   "Number_exons","Exons_Starts","Exons_Ends","NA2","Gene_name",
                                   "NA3","NA4","NA5"))

refseqbed <- read.table(refseqbed_path,h=F,sep="\t",stringsAsFactors=FALSE)
refseqshort <- read.table(refseqshort_path,h=F,sep="\t",stringsAsFactors=FALSE)

TTseq_list = list.files(TTseq_path, pattern="*.bed")
TTseq_list

for (i in seq_along(TTseq_list)) {
  filename <- sub(".bed", "", TTseq_list[i])
  filename <- sub("-", "_", filename)
  df <- read.table(paste0(TTseq_path,TTseq_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",paste0("TTseq_",filename)))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(filename, df)
}
nrow(MG9_11)
MG9_11[duplicated(MG9_11$Gene_name),]

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
MG9_11_reads <- 143598330
MG9_12_reads <- 150557906
MG9_13_reads <- 73388259
MG9_14_reads <- 131822641
MG9_15_reads <- 136900441
MG9_16_reads <- 124134871
MG9_17_reads <- 79232790
MG9_18_reads <- 85923486

# ------------------
# 1. Normalization |
# ------------------
MG9_11$TTseq_MG9_11 = (MG9_11$TTseq_MG9_11/MG9_11_reads) / (MG9_11$End-MG9_11$Start) * 10000000
MG9_12$TTseq_MG9_12 = (MG9_12$TTseq_MG9_12/MG9_12_reads) / (MG9_12$End-MG9_12$Start) * 10000000
MG9_13$TTseq_MG9_13 = (MG9_13$TTseq_MG9_13/MG9_13_reads) / (MG9_13$End-MG9_13$Start) * 10000000
MG9_14$TTseq_MG9_14 = (MG9_14$TTseq_MG9_14/MG9_14_reads) / (MG9_14$End-MG9_14$Start) * 10000000
MG9_15$TTseq_MG9_15 = (MG9_15$TTseq_MG9_15/MG9_15_reads) / (MG9_15$End-MG9_15$Start) * 10000000
MG9_16$TTseq_MG9_16 = (MG9_16$TTseq_MG9_16/MG9_16_reads) / (MG9_16$End-MG9_16$Start) * 10000000
MG9_17$TTseq_MG9_17 = (MG9_17$TTseq_MG9_17/MG9_17_reads) / (MG9_17$End-MG9_17$Start) * 10000000
MG9_18$TTseq_MG9_18 = (MG9_18$TTseq_MG9_18/MG9_18_reads) / (MG9_18$End-MG9_18$Start) * 10000000

sum(MG9_11$TTseq_MG9_11)
sum(MG9_13$TTseq_MG9_13)

sum(MG9_12$TTseq_MG9_12)
sum(MG9_14$TTseq_MG9_14)

sum(MG9_15$TTseq_MG9_15)
sum(MG9_17$TTseq_MG9_17)

sum(MG9_16$TTseq_MG9_16)
sum(MG9_18$TTseq_MG9_18)

# -----------
# Save data |
# -----------
write.table(MG9_11, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_11.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)     
write.table(MG9_12, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_12.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_13, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_13.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_14, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_14.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_15, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_15.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_16, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_16.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_17, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_17.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_18, 
            file = paste0(output,"Intersect_RefSeq_Normalized_MG9_18.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  

MG9_11$Gene_name=MG9_12$Gene_name
nrow(MG9_11)
summary(MG9_12)
summary(MG9_13)
summary(MG9_14)
summary(MG9_15)
summary(MG9_17)
# ---------------------------
# Determine TTseq threshold |
# ---------------------------

plot(density(MG9_15$TTseq_MG9_15, to=0.001),col="indianred1",lwd=2,main="reads_normalized")
lines(density(MG9_17$TTseq_MG9_17, to=0.001),lwd=2,col="indianred1")
lines(density(MG9_11$TTseq_MG9_11, to=0.001),lwd=2,col="deepskyblue")
lines(density(MG9_13$TTseq_MG9_13, to=0.001),lwd=2,col="deepskyblue")
legend("topright", c("WT1 0","WT2 0", "TKO1 0","TKO2 0"), col=c("deepskyblue", "deepskyblue", "indianred1","indianred1"), lwd=2, inset = .02)
abline(v=0.00006, col="grey2",lwd=1.5, lty=2)

plot(density(MG9_18$TTseq_MG9_18, to=0.005),col=2,lwd=2,main="reads_normalized")
lines(density(MG9_14$TTseq_MG9_14, to=0.005),lwd=2,col=4)
lines(density(MG9_16$TTseq_MG9_16, to=0.005),lwd=2,col=2)
lines(density(MG9_12$TTseq_MG9_12, to=0.005),lwd=2,col=4)
legend("topright", c("WT1 5","WT2 5", "TKO1 5","TKO2 5"), col=c("blue", "blue", "red1","red1"), lwd=2, inset = .02)
abline(v=0.0002, col="grey2",lwd=1.5, lty=2) # 0.002 as an optimal value to consider a threshold

# --------------
# 2. Filtering | Maslon doesn't say anything. Fuchs uses (600) as the maximum length of fragments aligned to consider
# --------------
MG9_11 <- subset(MG9_11, MG9_11$TTseq_MG9_11 > 0.002 & (MG9_11$End-MG9_11$Start)>20)
MG9_12 <- subset(MG9_12, MG9_12$TTseq_MG9_12 > 0.002 & (MG9_12$End-MG9_12$Start)>20)
MG9_13 <- subset(MG9_13, MG9_13$TTseq_MG9_13 > 0.002 & (MG9_13$End-MG9_13$Start)>20)
MG9_14 <- subset(MG9_14, MG9_14$TTseq_MG9_14 > 0.002 & (MG9_14$End-MG9_14$Start)>20)
MG9_15 <- subset(MG9_15, MG9_15$TTseq_MG9_15 > 0.002 & (MG9_15$End-MG9_15$Start)>20)
MG9_16 <- subset(MG9_16, MG9_16$TTseq_MG9_16 > 0.002 & (MG9_16$End-MG9_16$Start)>20)
MG9_17 <- subset(MG9_17, MG9_17$TTseq_MG9_17 > 0.002 & (MG9_17$End-MG9_17$Start)>20)
MG9_18 <- subset(MG9_18, MG9_18$TTseq_MG9_18 > 0.002 & (MG9_18$End-MG9_18$Start)>20)

# ----------------------
# Create common column |
# ----------------------
MG9_11$common <- paste0(MG9_11$Start, ";", MG9_11$End, ";", MG9_11$Gene_name)
MG9_12$common <- paste0(MG9_12$Start, ";", MG9_12$End, ";", MG9_12$Gene_name)
MG9_13$common <- paste0(MG9_13$Start, ";", MG9_13$End, ";", MG9_13$Gene_name)
MG9_14$common <- paste0(MG9_14$Start, ";", MG9_14$End, ";", MG9_14$Gene_name)
MG9_15$common <- paste0(MG9_15$Start, ";", MG9_15$End, ";", MG9_15$Gene_name)
MG9_16$common <- paste0(MG9_16$Start, ";", MG9_16$End, ";", MG9_16$Gene_name)
MG9_17$common <- paste0(MG9_17$Start, ";", MG9_17$End, ";", MG9_17$Gene_name)
MG9_18$common <- paste0(MG9_18$Start, ";", MG9_18$End, ";", MG9_18$Gene_name)

head(refseq)
refseq$common <- paste0(refseq$Start, ";", refseq$End, ";", refseq$Gene_name)

# ------------
# Merge data |
# ------------
#merged_TTseq_WT0<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F, all = T),
#                          list(MG9_11,
#                               MG9_13[,c("TTseq_MG9_13","common")]))

#merged_TTseq_WT5<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F, all = T),
#                          list(MG9_12,
#                               MG9_14[,c("TTseq_MG9_14","common")]))

#merged_TTseq_TKO0<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F, all = T),
#                          list(MG9_15,
#                               MG9_17[,c("TTseq_MG9_17","common")]))

#merged_TTseq_TKO5<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F, all = T),
#                          list(MG9_16,
#                               MG9_18[,c("TTseq_MG9_18","common")]))

merged_TTseq_all<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F, all = T),
                           list(MG9_11,
                                MG9_12[,c("TTseq_MG9_12","common")],
                                MG9_13[,c("TTseq_MG9_13","common")],
                                MG9_14[,c("TTseq_MG9_14","common")],
                                MG9_15[,c("TTseq_MG9_15","common")],
                                MG9_16[,c("TTseq_MG9_16","common")],
                                MG9_17[,c("TTseq_MG9_17","common")],
                                MG9_18[,c("TTseq_MG9_18","common")]))

merged_TTseq_all2<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F),
                          list(MG9_11,
                               MG9_12[,c("TTseq_MG9_12","common")],
                               MG9_13[,c("TTseq_MG9_13","common")],
                               MG9_14[,c("TTseq_MG9_14","common")],
                               MG9_15[,c("TTseq_MG9_15","common")],
                               MG9_16[,c("TTseq_MG9_16","common")],
                               MG9_17[,c("TTseq_MG9_17","common")],
                               MG9_18[,c("TTseq_MG9_18","common")]))
nrow(merged_TTseq_all)
tail(merged_TTseq_all2)
nrow(MG9_11)
nrow(MG9_13)
head(merged_TTseq_WT0, 20L)


dim(MG9_11[duplicated(MG9_11$common),])

merged_refseq<- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = F, all = T),
                       list(refseq,
                            merged_TTseq_all[,c("common")]))
# -----------
# Save data |
# -----------
write.table(merged_TTseq_all, 
            file = paste0(output,"Intersect_RefSeq_TTseq_all_Normalized.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)     
