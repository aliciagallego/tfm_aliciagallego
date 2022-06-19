#!/usr/bin/env Rscript

#############################
## Normalize (RefSeq 20621) # 20220329
#############################

# The objective of this script is the normalization by total number of reads per experiment of the transcriptome RefSeq list based on TTseq reads

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/2.2_Intersect_RefSeqBED3Kb_BAM/"
output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/3.2_Normalization_3Kb/"

# -----------
# Open data |
# -----------
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
nrow(MG9_18)
head(MG9_18)
MG9_12[duplicated(MG9_12$Gene_name),]

# ------------------------------------
# Total reads (data from experiment) | Only 5 min samples are taken into account
# ------------------------------------
#MG9_11_reads <- 143598330
#MG9_13_reads <- 73388259
#MG9_15_reads <- 136900441
#MG9_17_reads <- 79232790

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

sum(MG9_12$TTseq_MG9_12)
sum(MG9_14$TTseq_MG9_14)
sum(MG9_16$TTseq_MG9_16)
sum(MG9_18$TTseq_MG9_18)

head(MG9_12)
max(MG9_12$TTseq_MG9_12)

# --------------
# Size density | No aporta info relevante
# --------------
options(scipen=F)
plot(density((MG9_12$End-MG9_12$Start), to=100000),col="orange",lwd=2,main="transcript size distribution \n 20621 genes")
legend("topright", c("20Kb","30Kb", "40Kb","50Kb","60Kb"), col=c("gray47", "gray57", "gray67","gray77","gray87"), lwd=1.5, lty=2, inset = .02)
abline(v=20000, col="gray47",lwd=1.5, lty=2)
abline(v=30000, col="gray57",lwd=1.5, lty=2)
abline(v=40000, col="gray67",lwd=1.5, lty=2)
abline(v=50000, col="gray77",lwd=1.5, lty=2)
abline(v=60000, col="gray87",lwd=1.5, lty=2)

# ------------------------
# Explore #reads density |
# ------------------------
plot(density(MG9_18$TTseq_MG9_18, to=20),col="indianred1",lwd=2,main="reads_normalized")
lines(density(MG9_16$TTseq_MG9_16, to=20),lwd=2,col="red")
lines(density(MG9_12$TTseq_MG9_12, to=20),lwd=2,col="blue")
lines(density(MG9_14$TTseq_MG9_14, to=20),lwd=2,col="deepskyblue")
legend("topright", c("WT1 5 (MG9_12)","WT2 5 (MG9_14)", "TKO1 5 (MG9_16)","TKO2 5 (MG9_18)","threshold at 4"), col=c("blue", "deepskyblue", "red","indianred1","grey"), 
       lwd=c(2,2,2,2,1.5), lty=c(1,1,1,1,2), inset = .02)
abline(v=4, col="grey",lwd=1.5, lty=2)

# --------------------------------------------------------------
# Exploring expression levels for finding an optimal threshold | Select genes with > 0.7 reads
# --------------------------------------------------------------
MG9_12_fil <- MG9_12[MG9_12$TTseq_MG9_12>0.7,]
nrow(MG9_12_fil)
summary(MG9_12)

MG9_14_fil <- MG9_14[MG9_14$TTseq_MG9_14<0.7,]
nrow(MG9_14_fil)
summary(MG9_14)

MG9_16_fil <- MG9_16[MG9_16$TTseq_MG9_16<0.7,]
nrow(MG9_16_fil)
summary(MG9_16)

MG9_18_fil <- MG9_18[MG9_18$TTseq_MG9_18<0.7,]
nrow(MG9_18_fil)
summary(MG9_18)

MG9_12_by_rate <- MG9_12[order(MG9_12$TTseq_MG9_12, decreasing = TRUE),] # sort by size
head(MG9_12_by_rate)
MG9_12[which(MG9_12$Gene_name == "Dpy30"), ]
MG9_14[which(MG9_14$Gene_name == "Dpy30"), ]
MG9_16[which(MG9_16$Gene_name == "Dpy30"), ]
MG9_18[which(MG9_18$Gene_name == "Dpy30"), ]

# -----------------------
# Plot elongation rates |
# -----------------------
#png(file = paste0(output, "5_Plots/Plot_rates_and_reads_WT.png"))
MG9_12_ordered <- MG9_12[order(MG9_12$TTseq_MG9_12,decreasing = F),] 
MG9_14_ordered <- MG9_14[order(MG9_14$TTseq_MG9_14,decreasing = F),] 
MG9_16_ordered <- MG9_16[order(MG9_16$TTseq_MG9_16,decreasing = F),] 
MG9_18_ordered <- MG9_18[order(MG9_18$TTseq_MG9_18,decreasing = F),] 

plot(log(1+MG9_12_ordered$TTseq_MG9_12),
     main="Elongation rates",
     xlab="# genes",
     ylab="elongation rates",
     lwd = 2,
     type = "l",
     col="blue")
lines(log(1+MG9_14_ordered$TTseq_MG9_14),
      col="lightblue",
      lwd = 2,
      type = "l")
lines(log(1+MG9_16_ordered$TTseq_MG9_16),
      col="red",
      lwd = 2,
      type = "l")
lines(log(1+MG9_18_ordered$TTseq_MG9_18),
      col="salmon",
      lwd = 2,
      type = "l")
legend("topleft", inset=.02, legend=c("Number of reads/1000","Elongation rates Kb/min"),
       col=c("lightblue","blue"), lty=c(NA, 1),pch=c(1, NA), lwd=c(1,2), cex=0.8)
dev.off()	

head(MG9_18)

# -----------
# Save data |
# -----------
write.table(MG9_12, 
            file = paste0(output,"Intersect_RefSeq3Kb_Normalized_MG9_12.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_14, 
            file = paste0(output,"Intersect_RefSeq3Kb_Normalized_MG9_14.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_16, 
            file = paste0(output,"Intersect_RefSeq3Kb_Normalized_MG9_16.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  
write.table(MG9_18, 
            file = paste0(output,"Intersect_RefSeq3Kb_Normalized_MG9_18.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)  