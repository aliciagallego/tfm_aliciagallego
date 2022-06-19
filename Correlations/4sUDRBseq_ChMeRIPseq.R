#!/usr/bin/env Rscript

#############################################################
## Correlations elongation rates (TTseq) and m6A (meRIPseq) #
#############################################################

# -----------
# Libraries |
# -----------
#library("ggpubr")

# -------
# Paths |
# -------
m6AWT_path <- "/media/cc/A/Alicia/NGS/meRIP/meRIP_output/Intersect_RefSeq_BAM_repeated1/3_Normalized_IPInput/WU_normTotalReads_IPInput.bed"
m6ATKO_path <- "/media/cc/A/Alicia/NGS/meRIP/meRIP_output/Intersect_RefSeq_BAM_repeated1/3_Normalized_IPInput/TU_normTotalReads_IPInput.bed"
rates_path <- "/media/cc/A/Alicia/NGS/RNApolII/RNApolII_output/Intersect_RefSeq_BAM_repeated1/5_Merged_with_TTseq/Inputs/TTseq-all_Elongation_RerSef_5min_merged_withSize_clean_XY.txt"
output <- "/media/cc/A/Alicia/NGS/Correlations/TTseq_m6A/"

# -----------
# Open data |
# -----------
m6AWT <- read.table(m6AWT_path,h=F,sep="\t",stringsAsFactors=FALSE,
                    col.names = c("Chr","Start","End","All_inf","NA","Strand","Gene_name","meRIP_WT"))
m6ATKO <- read.table(m6ATKO_path,h=F,sep="\t",stringsAsFactors=FALSE,
                     col.names = c("Chr","Start","End","All_inf","NA","Strand","Gene_name","meRIP_TKO"))
rates <- read.table(rates_path,h=T,sep="\t",stringsAsFactors=FALSE,
                    col.names = c("Gene_name","TTseq_WT","TTseq_TKO","RN","ID","Chr","Strand","Start","End",
                                  "CDS_start","CDS_end","Exon_number","Exon_starts","Exon_ends","Z","cmpl1",
                                  "cmpl2","SP","Size"))

# -------------
# Merge files |
# -------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(m6AWT[,c("Chr","Start","End","NA.","Strand","meRIP_WT","Gene_name")],
                     m6ATKO[,c("meRIP_TKO","Gene_name")],
                     rates[,c("ID","TTseq_WT","TTseq_TKO","Gene_name")]))

# ------------------
# Save merged file |
# ------------------
write.table(merged, file = paste0(output,"TTseq-meRIP_merged.bed"), 
            quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

# ---------------------------------------------------------------
# Remove outliers from meRIP data and log transform RNApol data |
# ---------------------------------------------------------------
Q <- quantile(merged$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(merged, merged$meRIP_WT > low & merged$meRIP_WT < up)

Q <- quantile(merged$meRIP_TKO, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_TKO)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_TKO <- subset(merged, merged$meRIP_TKO > low & merged$meRIP_TKO < up)

# ------------------
# Correlation test |
# ------------------
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson", "kendall", "spearman"))

# Spearman resulted to be more significant
png(file = paste0(output, "Spearman_TTseq_m6A_WT.png"))
ggscatter(saved_values_WT, x = "meRIP_WT", y = "TTseq_WT",
          color="lightblue",shape = 0,
          size=0.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "m6A levels", ylab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation m6A levels vs. elongation rates - WT (2314 genes)")
dev.off()

# Pearson resulted to be more significant
png(file = paste0(output, "Pearson_TTseq_m6A_TKO.png"))
ggscatter(saved_values_TKO, x = "meRIP_TKO", y = "TTseq_TKO",
          color="pink",shape = 0,
          size=0.5,
          add.params = list(color = "lightcoral", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "m6A levels", ylab = "Elongation rate (Kb/min)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Pearson correlation m6A levels vs. elongation rates - TKO (2284 genes)")
dev.off()
