#!/usr/bin/env Rscript

##############################################
## Multiple correlation m6A, RNAPII, cheRNAs #
##############################################

# -----------
# Libraries |
# -----------
library(GGally)
library(psych)
library(ggpairs)
library(colorspace)

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull.txt"
refseq_path <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/4_Input_genes/Input_genes_20Kb.txt")

m6AWT_path <- "/media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/3_Normalized_data/meRIP_RefSeqLongList_Normalized_WT_2Kb.txt"
m6ATKO_path <- "/media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/3_Normalized_data/meRIP_RefSeqLongList_Normalized_TKO_2Kb.txt"

RNApolGB_path <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Intersect_RefSeq_BAM/2_Normalized_data/GeneBody/RNAPII_RefSeq_LongList_GB_500bp.bed"
RNApolPR_path <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Intersect_RefSeq_BAM/2_Normalized_data/Promoters/RNAPII_RefSeq_LongList_PR_500bp.bed"

Pausing_path <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Pausing_index/RNAPII_RefSeq_LongList_PausingIndex_500bp.txt"

cheRNA_path <- "/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/2_Normalized_data/cheRNA_RefSeq_LongList_Normalized.bed"

output <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/6_Correlations/")

# -----------
# Open data |
# -----------
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=FALSE)
refseq <- read.table(refseq_path,h=T,sep="\t",stringsAsFactors=FALSE)

m6AWT <- read.table(m6AWT_path, h=T,sep="\t",stringsAsFactors=F)
m6ATKO <- read.table(m6ATKO_path, h=T,sep="\t",stringsAsFactors=F)

RNApolGB <- read.table(RNApolGB_path, h=T, sep="\t", stringsAsFactors=F)
RNApolPR <- read.table(RNApolPR_path, h=T, sep="\t", stringsAsFactors=F)

Pausing <- read.table(Pausing_path,h=T,sep="\t",stringsAsFactors=F)

cheRNA <- read.table(cheRNA_path,h=T,sep="\t",stringsAsFactors=F)

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
                     TTseq[,c("Gene_name","WT_Pull","TK0_Pull")]))

# ---------------------------------------------------------
# Get transcript size in Kb based on End and Start coords |
# ---------------------------------------------------------
merged$Size <- (merged$End - merged$Start)/1000

# -------------------
# Delete NaN values |
# -------------------
merged_noNan <- na.omit(merged)  

# ------------------------------------------------
# Remove genes with elongation rate < 0.5 kb/min |
# ------------------------------------------------
merged_non05 <- merged_noNan[which(merged_noNan$WT_Pull > 0.5 & merged_noNan$TK0_Pull > 0.5),]
colnames(merged_non05)[1] <- "Gene_name"

# ----------------------------
# Merge rates and other data |
# ----------------------------
merged_WT_alldata<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                           list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                                RNApolGB[,c("RNAPII_WTmean_GB", "Gene_name")],
                                RNApolPR[,c("RNAPII_WTmean_PR", "Gene_name")],
                                Pausing[,c("PausIndex_WT","Gene_name")],
                                cheRNA[,c("cheRNA_WT_mean","Gene_name")],
                                merged_non05[,c("WT_Pull","Gene_name")]))

merged_TKO_alldata<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                            list(m6ATKO[,c("Chr","Start","End","Strand","meRIP_TKO","Gene_name")],
                                 RNApolGB[,c("RNAPII_TKOmean_GB", "Gene_name")],
                                 RNApolPR[,c("RNAPII_TKOmean_PR", "Gene_name")],
                                 Pausing[,c("PausIndex_TKO","Gene_name")],
                                 cheRNA[,c("cheRNA_TKO_mean","Gene_name")],
                                 merged_non05[,c("TK0_Pull","Gene_name")]))
# -----------------
# Remove outliers |
# -----------------

# meRIP
Q <- quantile(merged_WT_alldata$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_WT_alldata$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_WT_alldata <- subset(merged_WT_alldata, merged_WT_alldata$meRIP_WT > low & merged_WT_alldata$meRIP_WT < up)
nrow(merged_WT_alldata)

Q <- quantile(merged_TKO_alldata$meRIP_TKO, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_TKO_alldata$meRIP_TKO)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_TKO_alldata <- subset(merged_TKO_alldata, merged_TKO_alldata$meRIP_TKO > low & merged_TKO_alldata$meRIP_TKO < up)
nrow(merged_TKO_alldata)

# cheRNA
Q <- quantile(merged_WT_alldata$cheRNA_WT_mean, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_WT_alldata$cheRNA_WT_mean)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_WT_alldata <- subset(merged_WT_alldata, merged_WT_alldata$cheRNA_WT_mean > low & merged_WT_alldata$cheRNA_WT_mean < up)
nrow(merged_WT_alldata)

Q <- quantile(merged_TKO_alldata$cheRNA_TKO_mean, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_TKO_alldata$cheRNA_TKO_mean)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_TKO_alldata <- subset(merged_TKO_alldata, merged_TKO_alldata$cheRNA_TKO_mean > low & merged_TKO_alldata$cheRNA_TKO_mean < up)
nrow(merged_TKO_alldata)

# RNAPII_GB
Q <- quantile(merged_WT_alldata$RNAPII_WTmean_GB, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_WT_alldata$RNAPII_WTmean_GB)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_WT_alldata <- subset(merged_WT_alldata, merged_WT_alldata$RNAPII_WTmean_GB > low & merged_WT_alldata$RNAPII_WTmean_GB < up)
nrow(merged_WT_alldata)

Q <- quantile(merged_TKO_alldata$RNAPII_TKOmean_GB, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_TKO_alldata$RNAPII_TKOmean_GB)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_TKO_alldata <- subset(merged_TKO_alldata, merged_TKO_alldata$RNAPII_TKOmean_GB > low & merged_TKO_alldata$RNAPII_TKOmean_GB < up)
nrow(merged_TKO_alldata)

# RNAPII_PR
Q <- quantile(merged_WT_alldata$RNAPII_WTmean_PR, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_WT_alldata$RNAPII_WTmean_PR)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_WT_alldata <- subset(merged_WT_alldata, merged_WT_alldata$RNAPII_WTmean_PR > low & merged_WT_alldata$RNAPII_WTmean_PR < up)
nrow(merged_WT_alldata)

Q <- quantile(merged_TKO_alldata$RNAPII_TKOmean_PR, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged_TKO_alldata$RNAPII_TKOmean_PR)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
merged_TKO_alldata <- subset(merged_TKO_alldata, merged_TKO_alldata$RNAPII_TKOmean_PR > low & merged_TKO_alldata$RNAPII_TKOmean_PR < up)
nrow(merged_TKO_alldata)

# --------------
# Pairs panels |
# --------------

x <- matrix(rnorm(120*5),ncol=5)
col <- (rainbow_hcl(2))[c(rep(1, 60), rep(2,60))]
pairs(merged_WT_alldata[,c(6,7,8,10)], col = col, lower.panel = NULL, cex.labels=2, pch=19, cex = 0.8, cex.axis = 2)

colnames(merged_WT_alldata)[c(6,7,8,10)] <- c("m6A","RNAPII PR","RNAPII GB", "cheRNA levels")
colnames(merged_TKO_alldata)[c(6,7,8,10)] <- c("m6A","RNAPII PR","RNAPII GB", "cheRNA levels")

cor.test(merged_WT_alldata$m6A, merged_WT_alldata$RNAPII_PR, method="spearman")$p.val
cor.test(merged_WT_alldata$m6A, merged_WT_alldata$RNAPII_GB, method="spearman")$p.val
cor.test(merged_WT_alldata$m6A, merged_WT_alldata$cheRNA_levels, method="spearman")$p.val
cor.test(merged_WT_alldata$RNAPII_PR, merged_WT_alldata$cheRNA_levels, method="spearman")$p.val

cor.test(merged_TKO_alldata$m6A, merged_TKO_alldata$RNAPII_PR, method="spearman")$p.val
cor.test(merged_TKO_alldata$m6A, merged_TKO_alldata$RNAPII_GB, method="spearman")$p.val
cor.test(merged_TKO_alldata$m6A, merged_TKO_alldata$cheRNA_levels, method="spearman")$p.val
cor.test(merged_TKO_alldata$RNAPII_PR, merged_TKO_alldata$cheRNA_levels, method="spearman")$p.val

nrow(merged_WT_alldata)#4573
png(file = paste0(output, "/Correlation_panel_Pull_WT.png"))
pairs.panels(merged_WT_alldata[c(6,7,8,10)],
             density=T,
             pch=20,
             smoother=F,
             lm=T,
             method="spearman",
             hist.col=4,
             breaks = 20,
             rug=F,
             stars=T,
             cex.labels=1.5,
             ci=T,
             digits=3, cex=1)
dev.off()

nrow(merged_TKO_alldata)#4517
png(file = paste0(output, "/Correlation_panel_Pull_TKO.png"))
pairs.panels(merged_TKO_alldata[c(6,7,8,10)],
             density=T,
             smoother=F,
             col="red",
             pch=20,
             lm=T,
             method="spearman",
             hist.col=2,
             breaks = 20,
             rug=F,
             stars=T,
             cex.labels=1.5,
             ci=T,
             digits=3, cex=1)
dev.off()


ggpairs(merged_TKO_alldata, columns = c(6,7,8,10), lower=list(continuous="smooth"),
        diag = list(continuous = wrap("barDiag", colour = "blue")),
        upper = list(continuous = wrap("cor", size = 5))) 
