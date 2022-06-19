#!/usr/bin/env Rscript

########################################################
## Normalize H1 Cao et al 2013 data (RefSeq long list) # 20220609
########################################################

# -------
# Paths |
# -------
H1c_path <- "/media/cc/B/Alicia/H1_Cao/H1_Cao_output/1_bedtools_intersect/H1c_2kbTSS.bed"
H1d_path <- "/media/cc/B/Alicia/H1_Cao/H1_Cao_output/1_bedtools_intersect/H1d_2kbTSS.bed"
repli_path <- "/media/cc/B/Alicia/H1_Cao/H1_Cao_output/1_bedtools_intersect/SRR828743_2kbTSS.bed"
output <- "/media/cc/B/Alicia/H1_Cao/H1_Cao_output/"

# -----------
# Open data |
# -----------
H1c <- read.table(H1c_path,h=F,sep="\t",stringsAsFactors=F)
H1d <- read.table(H1d_path,h=F,sep="\t",stringsAsFactors=F)
repli <- read.table(repli_path,h=F,sep="\t",stringsAsFactors=F)
tail(H1d)

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
H1c_reads <-	55791364
H1d_reads <-	150519461

# ---------------
# Normalization |
# ---------------

H1c$V8 = H1c$V7/H1c_reads*10e06
H1d$V8 = H1d$V7/H1d_reads*10e06

# Colnames
names_vec <- c("Chr", "Start", "End", "Gene_name", "NA1","Strand", "readsCao","readsCaonorm")

colnames(H1c) <- names_vec
colnames(H1d) <- names_vec

head(H1c)

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull.txt"
refseq_path <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/4_Input_genes/Input_genes_20Kb.txt")

# -----------
# Open data |
# -----------
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=FALSE)
refseq <- read.table(refseq_path,h=T,sep="\t",stringsAsFactors=FALSE)

# ----------------------------------------------------
# Transform chr20, chr21 to chrX,chrY in refseq list | 
# ----------------------------------------------------
count <- 0
for(i in 1:nrow(refseq)) {
  if(refseq[i,3] == '20') {
    refseq[i,3]="X"
  }
  else if(refseq[i,3] == '21'){
    refseq[i,3]="Y"
  }
  count = count+1
}
print(count)

# -------------------------------------------------
# Remove white spaces in Gene_names in TTseq list | 
# -------------------------------------------------
TTseq$Gene_name <- gsub("\\s", "", TTseq$Gene_name)
TTseq[ TTseq == "NaN" ] <- NA

# ------------------------------
# Rates from bp/5min to Kb/min |
# ------------------------------
TTseq[,2] = TTseq[,2]/5000
TTseq[,3] = TTseq[,3]/5000

# ---------------------------------------------------
# Merge Elongation rates list and Input RefSeq list |
# ---------------------------------------------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by.x="Name", by.y="Gene_name", sort = F),
                list(refseq[,c("Chromosome","Start", "End", "Exon_number", "Orientation","Name")],
                     TTseq[,c("Gene_name","WT_Pull", "TK0_Pull")]))

nrow(merged)
nrow(TTseq)
nrow(refseq)
# ---------------------------------------------------------
# Get transcript size in Kb based on End and Start coords |
# ---------------------------------------------------------
merged$Size <- (merged$End - merged$Start)/1000

# -------------------
# Delete NaN values |
# -------------------
merged_noNan <- na.omit(merged)  
nrow(merged_noNan)

# ------------------------------------------------
# Remove genes with elongation rate < 0.5 kb/min |
# ------------------------------------------------
nrow(merged_noNan)
nrow(merged_non05)
merged_non05 <- merged_noNan[which(merged_noNan$WT_Pull > 0.5 & merged_noNan$TK0_Pull > 0.5),]
head(merged_noNan)
colnames(merged_non05)[1] <- "Gene_name"

# ----------------------------
# Merge rates and other data |
# ----------------------------
merged_H1c<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                           list(H1c[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                                merged_non05[,c("WT_Pull", "Exon_number", "Size","Gene_name")]))

merged_H1d<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                    list(H1d[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                         merged_non05[,c("WT_Pull", "Exon_number", "Size","Gene_name")]))

head(H1c)
nrow(H1d)
nrow(merged_H1c)
nrow(merged_H1d)
head(merged_non05)

# ------------------------------
# Obtain tertiles and classify |
# ------------------------------

# WT
# --
QWT <- quantile(merged_non05$WT_Pull, probs = seq(0, 1, 1/3), na.rm = FALSE)  
upWT <-  QWT[[3]]  
lowWT<- QWT[[2]]

WTslow <- (subset(merged_non05, merged_non05$WT_Pull<lowWT))
WTslow <- WTslow[order(WTslow$WT_Pull),]
WTmedium  <- (subset(merged_non05, merged_non05$WT_Pull>=lowWT & merged_non05$WT_Pull<=upWT))
WTmedium <- WTmedium[order(WTmedium$WT_Pull),]
WTfast <- (subset(merged_non05, merged_non05$WT_Pull>upWT))
WTfast <- WTfast[order(WTfast$WT_Pull),]

# TKO
# ---
QTKO <- quantile(merged_non05$TK0_Pull, probs = seq(0, 1, 1/3), na.rm = FALSE)  
upTKO <-  QTKO[[3]]  
lowTKO<- QTKO[[2]]

TKOslow <- (subset(merged_non05, merged_non05$TK0_Pull<lowTKO))
TKOslow <- TKOslow[order(TKOslow$TK0_Pull),]
TKOmedium  <- (subset(merged_non05, merged_non05$TK0_Pull>=lowTKO & merged_non05$TK0_Pull<=upTKO))
TKOmedium <- TKOmedium[order(TKOmedium$TK0_Pull),]
TKOfast <- (subset(merged_non05, merged_non05$TK0_Pull>upTKO))
TKOfast <- TKOfast[order(TKOfast$TK0_Pull),]

# --------------------------------
# Merged tertiles and other data |
# --------------------------------
names_vec <- c("Gene_name", "Chr", "Start", "End", "Exon_number","Strand", "TTseq_WT_pull", 
               "TTseq_TKO_pull","Size")

colnames(WTslow) <- names_vec
colnames(WTmedium) <- names_vec
colnames(WTfast) <- names_vec
colnames(TKOslow) <- names_vec
colnames(TKOmedium) <- names_vec
colnames(TKOfast) <- names_vec

# merged

merged_WTslow_H1c<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(H1c[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                                      WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTmed_H1c<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(H1c[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                                WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTfast_H1c<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                          list(H1c[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                               WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTslow_H1d<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(H1d[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                                WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTmed_H1d<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                          list(H1d[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                               WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))
merged_WTfast_H1d<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(H1d[,c("Chr","Start","End","Strand","readsCaonorm","Gene_name")],
                                WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

# ----------
# Boxplots |
# ----------

png(file = paste0(output, "BoxPlots_H1Cao_tertiles.png"))
wilcox <- wilcox.test(merged_WTslow_H1c$readsCaonorm, merged_WTmed_H1c$readsCaonorm)
wilcox
wilcox <- wilcox.test(merged_WTmed_H1c$readsCaonorm,merged_WTfast_H1c$readsCaonorm)
wilcox
wilcox <- wilcox.test(merged_WTslow_H1c$readsCaonorm, merged_WTfast_H1c$readsCaonorm)
wilcox
boxplot(merged_WTslow_H1c$readsCaonorm, merged_WTmed_H1c$readsCaonorm,merged_WTfast_H1c$readsCaonorm,
        merged_WTslow_H1d$readsCaonorm, merged_WTmed_H1d$readsCaonorm,merged_WTfast_H1d$readsCaonorm,
        col = "grey",
        at = c(1,2,3,5,6,7),
        names = c("H1c_slow","H1c_med", "H1c_fast", "H1c_slow","H1c_med", "H1c_fast"),
        main = "H1 data from Cao et al.",
        ylab="H1 at TSS",
        outline = F,
        las=2,
        ylim=c(1,35),
        boxwex = 0.7)
segments(1,28,2,28,lwd=0.8)
segments(1,31,3,31,lwd=0.8)
segments(5,28,6,28,lwd=0.8)
segments(5,31,7,31,lwd=0.8)

text(1.5, 29,"*",cex=2,col=("black"))
text(2, 32, "*",cex=2,col=("black"))
text(5.5, 29,"*",cex=2,col=("black"))
text(6, 32, "*",cex=2,col=("black"))
dev.off()	
