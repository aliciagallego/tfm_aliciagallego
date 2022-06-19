#!/usr/bin/env Rscript

#####################################################
## Correlations RNApolII and m6A (RefSeq Long List) # 20220214
#####################################################

# -----------
# Libraries |
# -----------
library("ggpubr")

# -------
# Paths |
# -------
RNApolGB_path <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Intersect_RefSeq_BAM/2_Normalized_data/GeneBody/RNAPII_RefSeq_LongList_GB_500bp.bed"
RNApolPR_path <- "/media/cc/A/Alicia/NGS/RNApolII_3/RNApolII3_output/Intersect_RefSeq_BAM/2_Normalized_data/Promoters/RNAPII_RefSeq_LongList_PR_500bp.bed"

m6AWT_path <- "/media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/3_Normalized_data/meRIP_RefSeqLongList_Normalized_WT_2Kb.txt"
m6ATKO_path <- "/media/cc/A/Alicia/NGS/meRIP_2/meRIP2_output/3_Normalized_data/meRIP_RefSeqLongList_Normalized_TKO_2Kb.txt"

output <- "/media/cc/A/Alicia/NGS/Correlations/RNApol_m6A_RefSeqLongList/"

# -----------
# Open data |
# -----------
RNApolGB <- read.table(RNApolGB_path, h=T, sep="\t", stringsAsFactors=F)
RNApolPR <- read.table(RNApolPR_path, h=T, sep="\t", stringsAsFactors=F)

m6AWT <- read.table(m6AWT_path, h=T,sep="\t",stringsAsFactors=F)
m6ATKO <- read.table(m6ATKO_path, h=T,sep="\t",stringsAsFactors=F)

head(RNApolGB)
head(RNApolPR)

# -------------
# Merge files |
# -------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                     m6ATKO[,c("meRIP_TKO","Gene_name")],
                     RNApolGB[,c("RNAPII_WTmean_GB","RNAPII_TKOmean_GB", "Gene_name")],
                     RNApolPR[,c("RNAPII_WTmean_PR","RNAPII_TKOmean_PR", "Gene_name")]))

colnames(merged) <- c("Gene_name","Chr","Start","End","Strand","meRIP_WT","meRIP_TKO",
                      "RNApol_WT_GB","RNApol_TKO_GB","RNApol_WT_PR","RNApol_TKO_PR")

nrow(merged)
# ------------------
# Save merged file |
# ------------------
write.table(merged, file = paste0(output,"RNApol-meRIP_merged_RefSeqLongList_500bp.bed"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# ---------------------------------------------------------------
# Remove outliers from meRIP data and log transform RNApol data |
# ---------------------------------------------------------------
Q <- quantile(merged$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(merged, merged$meRIP_WT > low & merged$meRIP_WT < up)
head(saved_values_WT)
nrow(saved_values_WT)

saved_values_WT$RNApol_WT_GB <- log(saved_values_WT$RNApol_WT_GB*10000000+1)
saved_values_WT$RNApol_WT_PR <- log(saved_values_WT$RNApol_WT_PR*10000000+1)

Q <- quantile(merged$meRIP_TKO, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_TKO)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_TKO <- subset(merged, merged$meRIP_TKO > low & merged$meRIP_TKO < up)
head(saved_values_TKO)
nrow(saved_values_TKO)

saved_values_TKO$RNApol_TKO_GB <- log(saved_values_TKO$RNApol_TKO_GB*10000000+1)
saved_values_TKO$RNApol_TKO_PR <- log(saved_values_TKO$RNApol_TKO_PR*10000000+1)

# ------------------
# Correlation test |
# ------------------
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson", "kendall", "spearman"))
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson"))

# Spearman resulted to be more significant
png(file = paste0(output, "Spearman_RNApolII_meRIP_WT_GB_500bp.png"))
ggscatter(saved_values_WT, y = "meRIP_WT", x = "RNApol_WT_GB",
          color="lightblue",shape = 17,
          size=1.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "m6A levels", xlab = "log RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation RNApolII vs. m6A levels \n WT GB (15939 genes, RefSeq Long List)")
dev.off()

png(file = paste0(output, "Spearman_RNApolII_meRIP_TKO_GB_GB_500bp.png"))
ggscatter(saved_values_TKO, y = "meRIP_TKO", x = "RNApol_TKO_GB",
          color="pink",shape = 17,
          size=1.5,
          add.params = list(color = "lightcoral", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "m6A levels", xlab = "log RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation RNApolII vs. m6A levels \n TKO GB (15968 genes, RefSeq Long List)")
dev.off()

png(file = paste0(output, "Spearman_RNApolII_meRIP_WT_PR_GB_500bp.png"))
ggscatter(saved_values_WT, y = "meRIP_WT", x = "RNApol_WT_PR",
          color="steelblue1",shape = 17,
          size=1.5,
          add.params = list(color = "steelblue3", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "m6A levels", xlab = "log RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation RNApolII vs. m6A levels \n WT PR (15939 genes, RefSeq Long List)")
dev.off()

png(file = paste0(output, "Spearman_RNApolII_meRIP_TKO_PR_GB_500bp.png"))
ggscatter(saved_values_TKO, y = "meRIP_TKO", x = "RNApol_TKO_PR",
          color="red",shape = 17,
          size=1.5,
          add.params = list(color = "red3", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "m6A levels", xlab = "log RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation RNApolII vs. m6A levels \n TKO PR (15968 genes, RefSeq Long List)")
dev.off()