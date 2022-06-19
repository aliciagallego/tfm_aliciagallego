#!/usr/bin/env Rscript

#####################################################################
## Correlations elongation rates (4sUDRBseq or TTseq) and CheRNAseq #
#####################################################################

# -----------
# Libraries |
# -----------
library("ggpubr")
library("reshape2")

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull_processed_without05.txt"
cheRNA_path <- "/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/2_Normalized_data/cheRNA_RefSeq_LongList_Normalized.bed"

output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/6_Correlations/"

# -----------
# Open data |
# -----------
cheRNA <- read.table(cheRNA_path,h=T,sep="\t",stringsAsFactors=F)
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))

# -------------
# Merge files |
# -------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_WT_mean","cheRNA_TKO_mean","Gene_name")],
                     TTseq[,c("TTseq_WT_pull", "TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# ------------------
# Save merged file |
# ------------------
write.table(merged, file = paste0(output,"TTseq_cheRNA/TTseq_cheRNA.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# ----------------------------------
# Remove outliers from TTseqA data |
# ----------------------------------
Q <- quantile(merged$TTseqA_WT_TPM, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$TTseqA_WT_TPM)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(merged, merged$TTseqA_WT_TPM > low & merged$TTseqA_WT_TPM < up)
head(saved_values_WT)
nrow(saved_values_WT)

Q <- quantile(merged$TTseqA_TKO_TPM, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$TTseqA_TKO_TPM)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_TKO <- subset(merged, merged$TTseqA_TKO_TPM > low & merged$TTseqA_TKO_TPM < up)
head(saved_values_TKO)
nrow(saved_values_TKO)

# select TMP > 1
saved_values_WT2 <- subset(saved_values_WT, saved_values_WT$TTseqA_WT_TPM > 1)
saved_values_TKO2 <- subset(saved_values_TKO, saved_values_TKO$TTseqA_TKO_TPM > 1)

# select especific TMP based on TPM distributions
saved_values_WT3 <- subset(saved_values_WT, saved_values_WT$TTseqA_WT_TPM > 30)
saved_values_TKO3 <- subset(saved_values_TKO, saved_values_TKO$TTseqA_TKO_TPM > 15)

# ------------------
# Correlation test |
# ------------------
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson", "kendall", "spearman"))

# WT all reads
png(file = paste0(output, "TTseq_cheRNA/Spearman_TTseqPull_cheRNAs_WT_allreads.png"))
number <- nrow(merged)
ggscatter(merged, y = "cheRNA_WT_mean", x = "TTseq_WT_pull",
          color="blue",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "blue", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = T, cor.method = "spearman",
          ylim=c(0,1e-05),
          ylab = "Expression levels (cheRNA reads)", xlab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "blue", label.x = 6, label.y = 7.5e-06,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n expression levels (all reads) vs. elongation rates in WT \n (", number, " genes)"))
dev.off()

# TKO all reads
png(file = paste0(output, "TTseq_cheRNA/Spearman_TTseqPull_cheRNAs_TKO_allreads.png"))
number <- nrow(merged)
ggscatter(merged, x = "TTseq_TKO_pull", y = "cheRNA_TKO_mean",
          color="red2",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "red2", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = T, cor.method = "spearman",
          ylim=c(0,1e-05),
          ylab = "Expression levels (cheRNA reads)", xlab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "red2", label.x = 6, label.y = 7.5e-06,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n expression levels (all reads) vs. elongation rates in TKO \n (", number, " genes)"))
dev.off()