#!/usr/bin/env Rscript

#######################################################
## Correlations RNApolII and elongation rates (TTseq) # 20220119
#######################################################

# -----------
# Libraries |
# -----------
library("ggpubr")

# -------
# Paths |
# -------
RNApol_path <- "/media/cc/A/Alicia/NGS/RNApolII/RNApolII_output/Intersect_RefSeq_BAM_repeated1/5_Merged_with_TTseq/Inputs/"
output <- "/media/cc/A/Alicia/NGS/Correlations/RNApol_TTseq/"

# --------------------------------------------------------
# Read all files from directories and put them in a list |
# --------------------------------------------------------
GB_list = list.files(RNApol_path, pattern="*GB")
Prom_list = list.files(RNApol_path,pattern="promoter")
TT_list = list.files(RNApol_path, pattern="*withSize_clean_XY.txt")

# ----------------------------------------------------------------------------------
# Open RNApol and TTseq BED files, give colnames, create 'common' column and merge |
# ----------------------------------------------------------------------------------

# Gene Body
# ---------
for (i in seq_along(GB_list)) {
  filenameGB <- sub("_.*", "", GB_list[[i]])
  dfGB <- read.table(GB_list[[i]], header = FALSE, sep="\t", stringsAsFactors=FALSE, 
                     col.names = c("Chr","Start","End", "All_info", "NA", "Strand", "Gene_name","RNApol"))
  dfGB$common <- paste0(dfGB$Start, ";", dfGB$End, ";", dfGB$Gene_name)
  
  for (j in seq_along(TT_list)) {
    filenameTT <- sub("_.*", "", TT_list[[j]])
    dfTT <- read.table(TT_list[[j]], header = T, sep="\t", stringsAsFactors=FALSE, 
                       col.names = c("Gene_name","TTseq_WT","TTseq_TKO","RN","ID","Chr","Strand","Start","End",
                                     "CDS_start","CDS_end","Exon_number","Exon_starts","Exon_ends","Z","cmpl1",
                                     "cmpl2","SP","Size"))
    dfTT$common <- paste0(dfTT$Start, ";", dfTT$End, ";", dfTT$Gene_name)
    
    merged <- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = FALSE),
                     list(dfGB,
                          dfTT[,c("TTseq_WT","TTseq_TKO","common")])) # NOTE: I am keeping rates for WT and TKO
    merged <- merged[-c(1)]
    write.table(merged, file = paste0(output,filenameGB,"_",filenameTT,"_merged.bed"), 
                quote = F, sep="\t", col.names = TRUE, row.names = FALSE)
  }  
}

# Promoters
# ---------
for (i in seq_along(Prom_list)) {
  filenamePR <- sub("_.*", "", Prom_list[[i]])
  dfPR <- read.table(Prom_list[[i]], header = FALSE, sep="\t", stringsAsFactors=FALSE, 
                     col.names = c("Chr","Start","End", "All_info", "NA", "Strand", "Gene_name","RNApol"))
  dfPR$common <- paste0(dfPR$Gene_name)
  
  for (j in seq_along(TT_list)) {
    filenameTT <- sub("_.*", "", TT_list[[j]])
    dfTT <- read.table(TT_list[[j]], header = T, sep="\t", stringsAsFactors=FALSE, 
                       col.names = c("Gene_name","TTseq_WT","TTseq_TKO","RN","ID","Chr","Strand","Start","End",
                                     "CDS_start","CDS_end","Exon_number","Exon_starts","Exon_ends","Z","cmpl1",
                                     "cmpl2","SP","Size"))
    dfTT$common <- paste0(dfTT$Gene_name)
    
    merged <- Reduce(function(x,y) merge(x = x, y = y, by = "common", sort = FALSE),
                     list(dfPR,
                          dfTT[,c("TTseq_WT","TTseq_TKO","common")])) # NOTE: I am keeping rates for WT and TKO
    merged <- merged[-c(1)]
    write.table(merged, file = paste0(output,filenamePR,"_",filenameTT,"_merged.bed"), 
                quote = F, sep="\t", col.names = TRUE, row.names = FALSE)
    
  }  
}

# --------------------
# Open files - means |
# --------------------
WT_GB <- as.data.frame(read.table(paste0(output,"RNApol-GB-WTmean_TTseq-all_merged.bed"),
                                  header=T, sep="\t", stringsAsFactors=FALSE,
                                  col.names = c("Chr","Start","End","All_info","NA","Strand","Gene_name","RNApol_WT_GB",
                                                "TTseq_WT","TTseq_TKO")))

TKO_GB <- as.data.frame(read.table(paste0(output,"RNApol-GB-TKOmean_TTseq-all_merged.bed"),
                                                header=T, sep="\t", stringsAsFactors=FALSE,
                                   col.names = c("Chr","Start","End","All_info","NA","Strand","Gene_name","RNApol_TKO_GB",
                                                 "TTseq_WT","TTseq_TKO")))

WT_PR <- as.data.frame(read.table(paste0(output,"RNApol-Prom-WTmean_TTseq-all_merged.bed"),
                                  header=T, sep="\t", stringsAsFactors=FALSE,
                                  col.names = c("Chr","Start","End","All_info","NA","Strand","Gene_name","RNApol_WT_PR",
                                                "TTseq_WT","TTseq_TKO")))

TKO_PR <- as.data.frame(read.table(paste0(output,"RNApol-Prom-TKOmean_TTseq-all_merged.bed"),
                                   header=T, sep="\t", stringsAsFactors=FALSE,
                                   col.names = c("Chr","Start","End","All_info","NA","Strand","Gene_name","RNApol_TKO_PR",
                                                 "TTseq_WT","TTseq_TKO")))

# -------------
# Merge files |
# -------------
dim(TKO_GB[duplicated(WT_PR$Gene_name),]) # to confirm there are no duplications in gene names and we can merge four dataframes

merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                list(WT_GB[,c("Chr","Start","End","All_info","NA.","Strand","TTseq_WT","TTseq_TKO","RNApol_WT_GB","Gene_name")],
                                     TKO_GB[,c("RNApol_TKO_GB","Gene_name")],
                                     WT_PR[,c("RNApol_WT_PR", "Gene_name")],
                                     TKO_PR[,c("RNApol_TKO_PR", "Gene_name")]))
nrow(merged)
# ------------------
# Save merged file |
# ------------------
write.table(merged, file = paste0(output,"RNApol-TTseq_merged.bed"), 
            quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

# ------------------
# Correlation test |
# ------------------
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson", "kendall", "spearman"))
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson"))

# pearson resulted more significant for gene body
png(file = paste0(output, "Pearson_RNApolII_TTseq_WT_GB.png"))
ggscatter(merged, x = "TTseq_WT", y = "RNApol_WT_GB",
          color="lightblue",shape = 1,
          size=2,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Pearson correlation RNApolII vs. elongation rates - WT GB (2399 genes)")
dev.off()

png(file = paste0(output, "Pearson_RNApolII_TTseq_TKO_GB.png"))
ggscatter(merged, x = "TTseq_TKO", y = "RNApol_TKO_GB",
          color="pink",shape = 1,
          size=2,
          add.params = list(color = "lightcoral", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Pearson correlation RNApolII vs. elongation rates - TKO GB (2399 genes)")
dev.off()	

# spearman resulted more significant for promoters
png(file = paste0(output, "Spearman_RNApolII_TTseq_WT_PR.png"))
ggscatter(merged, x = "TTseq_WT", y = "RNApol_WT_PR",
          color="steelblue1",shape = 1,
          size=2,
          add.params = list(color = "steelblue3", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation RNApolII vs. elongation rates - WT PR (2399 genes)")
dev.off()	

png(file = paste0(output, "Spearman_RNApolII_TTseq_TKO_PR.png"))
ggscatter(merged, x = "TTseq_TKO", y = "RNApol_TKO_PR",
          color="red",shape = 1,
          size=2,
          add.params = list(color = "red3", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation RNApolII vs. elongation rates - TKO PR (2399 genes)")
dev.off()	

# -------------------------------
# Plot RNApol in the 2400 genes |
# -------------------------------
png(file = paste0(output, "Plot_RNApolII_2400_62_genes.png"))
merged_ordered <- merged[order(merged$RNApol_WT_GB,decreasing = F),] 
plot(merged_ordered$RNApol_WT_GB,
     col="lightblue",
     main="RNApolII occupancies (2300 genes)",
     xlab="transcripts",
     ylab="RNApolII",
     lwd = 2,
     type = "l")
merged_ordered <- merged[order(merged$RNApol_TKO_GB,decreasing = F),] 
lines(merged_ordered$RNApol_TKO_GB,
      col="pink",
      lwd = 1,
      type = "l")
merged_ordered <- merged[order(merged$RNApol_WT_PR,decreasing = F),] 
lines(merged_ordered$RNApol_WT_PR,
      col="blue",
      lwd = 1,
      type = "l")
merged_ordered <- merged[order(merged$RNApol_TKO_PR,decreasing = F),] 
lines(merged_ordered$RNApol_TKO_PR,
      col="red",
      lwd = 1,
      type = "l")
legend("topright", inset=.02, legend=c("WT  gene body", "WT  promotor","TKO gene body", "TKO promotor"),
       col=c("lightblue","blue", "pink", "red"), lty=1, lwd=2, cex=0.8)
dev.off()	