#!/usr/bin/env Rscript

###########################################################
## Normalize meRIP (Ensemble) and merge with Maslon rates #
###########################################################

# -----------
# Libraries |
# -----------
library("ggpubr")

# -------
# Paths |
# -------
m6A_path <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_output/2_Intersect_bedtools/"
rates_Maslon_path <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_data/Elongation_rates_MaslonWT.txt"
rates_out <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_output/"

# ---------------------------------------
# Open TXT data Maslon elongation rates |
# ---------------------------------------
rates_Maslon <- read.table(rates_Maslon_path,h=T,sep="\t",stringsAsFactors=FALSE)

# ---------------------------------------------
# Transform Maslon elongation rates to kb/min |
# ---------------------------------------------
rates_Maslon[,4] = rates_Maslon[,4]/1000
colnames(rates_Maslon) <- c("enst_id","R_squared","p_value","Maslon_Kb_min","tx_length","FPKM_control","boundary_5min",
                        "boundary_15min","gene_short_name")

# -------------------------------------
# Open BED data meRIP-Ensembl aligned |
# -------------------------------------
BED_list = list.files(m6A_path, pattern="*.bed")

for (i in seq_along(BED_list)) {
  filename <- sub("_2Kb.bed*", "", BED_list[i])
  df <- read.table(paste0(m6A_path,BED_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand",
                                 "NA1","NA2","Gene_name","Gene_type","meRIP"))
  df <- subset(df, select=-c(NA1,NA2)) # remove uninformative columns
  df <- df[!duplicated(df$Transcript.stable.ID),] # remove duplicated rows
  assign(filename, df)
}

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
WU1_IP_reads	<- 52603601
WU2_IP_reads	<- 56468821
TU1_IP_reads	<- 52209521
TU2_IP_reads	<- 53790404

WU1_In_reads	<- 35203047
WU2_In_reads	<- 36358279
TU1_In_reads	<- 34963985
TU2_In_reads	<- 34564725

# ------------------------------------
# Normalize by total reads and IP/In |
# ------------------------------------
WU1_IP$meRIP_norm = (WU1_IP$meRIP/WU1_IP_reads) / (WU1_In$meRIP/WU1_In_reads)
WU2_IP$meRIP_norm = (WU2_IP$meRIP/WU2_IP_reads) / (WU2_In$meRIP/WU2_In_reads)
TU1_IP$meRIP_norm = (TU1_IP$meRIP/TU1_IP_reads) / (TU1_In$meRIP/TU1_In_reads)
TU2_IP$meRIP_norm = (TU2_IP$meRIP/TU2_IP_reads) / (TU2_In$meRIP/TU2_In_reads)

# -----------------------------
# Remove NA, Inf and 0 values |
# -----------------------------
WU1_IP <- subset(WU1_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")
WU2_IP <- subset(WU2_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")
TU1_IP <- subset(TU1_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")
TU2_IP <- subset(TU2_IP, meRIP_norm!="NaN" & meRIP_norm!=0 & meRIP_norm!="Inf")

# -----------------
# Merge m6A files |
# -----------------

# WT
# --
mergedWT <- Reduce(function(x,y) merge(x = x, y = y, by = "Transcript.stable.ID", sort = F),
                   list(WU1_IP[,c("Chr","Start","End","Gene.stable.ID","Strand","Gene_name",
                                  "Gene_type","meRIP_norm","Transcript.stable.ID")],
                        WU2_IP[,c("meRIP_norm","Transcript.stable.ID")]))

mergedWT <- mergedWT[,c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand","Gene_name",
                        "Gene_type","meRIP_norm.x","meRIP_norm.y")]

colnames(mergedWT) <- c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand","Gene_name",
                        "Gene_type","meRIP_WT1","meRIP_WT2")

mergedWT$meRIP_WT <- (mergedWT$meRIP_WT1+mergedWT$meRIP_WT2)/2

# TKO
# ---
mergedTKO <- Reduce(function(x,y) merge(x = x, y = y, by = "Transcript.stable.ID", sort = F),
                   list(TU1_IP[,c("Chr","Start","End","Gene.stable.ID","Strand","Gene_name",
                                  "Gene_type","meRIP_norm","Transcript.stable.ID")],
                        TU2_IP[,c("meRIP_norm","Transcript.stable.ID")]))

mergedTKO <- mergedTKO[,c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand","Gene_name",
                        "Gene_type","meRIP_norm.x","meRIP_norm.y")]

colnames(mergedTKO) <- c("Chr","Start","End","Gene.stable.ID","Transcript.stable.ID","Strand","Gene_name",
                        "Gene_type","meRIP_TKO1","meRIP_TKO2")

mergedTKO$meRIP_TKO <- (mergedTKO$meRIP_TKO1+mergedTKO$meRIP_TKO2)/2

mergedWT[which(mergedWT$Transcript.stable.ID == "ENSMUST00000071543"), ]

# -----------
# Save data |
# -----------
write.table(mergedWT, 
            file = paste0(rates_out,"3_Normalized_data/meRIP_Ensembl_Normalized_WT_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

write.table(mergedTKO, 
            file = paste0(rates_out,"3_Normalized_data/meRIP_Ensembl_Normalized_TKO_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

# -------------------------------
# Merge m6A WT and Maslon files |
# -------------------------------
rates_Maslon_m6A_WT<- Reduce(function(x,y) merge(x = x, y = y, by.x = "Transcript.stable.ID",by.y = "enst_id", sort = F),
                             list(mergedWT[,c("Chr", "Start", "End", "Gene.stable.ID","Strand", "Gene_name", "Gene_type",
                                              "meRIP_WT1","meRIP_WT2","meRIP_WT","Transcript.stable.ID")],
                                  rates_Maslon[,c("Maslon_Kb_min","gene_short_name","enst_id")]))

lost <- nrow(rates_Maslon) - nrow(rates_Maslon_m6A_WT)
print(paste(lost,"transcripts lost from Maslon list after merge"))

# -----------
# Save data |
# -----------
write.table(rates_Maslon_m6A_WT, 
            file = paste0(rates_out,"4_meRIP_Maslon/meRIP_Normalized_Maslon_WT_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

# ---------------------
# Remove m6A outliers |
# ---------------------
Q <- quantile(rates_Maslon_m6A_WT$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(rates_Maslon_m6A_WT$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(rates_Maslon_m6A_WT, rates_Maslon_m6A_WT$meRIP_WT > low & rates_Maslon_m6A_WT$meRIP_WT < up)

lost2<-nrow(rates_Maslon)-nrow(saved_values_WT)
print(paste(lost2,"transcripts lost from Maslon list after merge and meRIP outlier removal"))

# -----------
# Save data |
# -----------
write.table(saved_values_WT, 
            file = paste0(rates_out,"4_meRIP_Maslon/meRIP_Normalized_Maslon_WT_noOutliers_2Kb.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

# ----------------------------
# Mann-Whitney-Wilcoxon test |
# ----------------------------
sink("meRIP_All_WilcoxTest.txt")
wilcox.test(mergedWT$meRIP_WT, mergedTKO$meRIP_TKO)
sink()

# ----------
# Boxplots |
# ----------
png(file = paste0(rates_out, "5_Plots/meRIP_All_2Kb"))
boxplot(mergedWT$meRIP_WT, mergedTKO$meRIP_TKO,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "m6A - means (Ensembl 119007 genes)",
        ylab="m6A levels",
        outline = F,
        xlab=paste("Wilcox test: W = 7866343962, p-value < 2.2e-16"),
        boxwex = 0.3)
dev.off()	

png(file = paste0(rates_out, "5_Plots/meRIP_replicates_2Kb"))
boxplot(mergedWT$meRIP_WT1, mergedWT$meRIP_WT2, mergedTKO$meRIP_TKO1,mergedTKO$meRIP_TKO2,
        col = c(4,4,2,2),
        at = c(1,2,4,5),
        names = c("WT1","WT2","TKO1", "TKO2"),
        main = "m6A - replicates (Ensembl)",
        ylab="m6A levels",
        outline = F,
        xlab=paste("(WT: 119007 genes)    (TKO: 119066 genes)"),
        boxwex = 0.5)
dev.off()

# -------------------------------
# Plots elongation rates Maslow |
# -------------------------------
rates_Maslon_m6A_WT <- rates_Maslon_m6A_WT[order(rates_Maslon_m6A_WT$Maslon_Kb_min,decreasing = TRUE),] 

png(file = paste0(rates_out, "5_Plots/Plot_MaslonElongation_meRIP_WT_2Kb.png"))
plot(rates_Maslon_m6A_WT$Maslon_Kb_min,
     col="blue",
     main="Elongation rates Maslon (808 genes + 2Kb)",
     xlab="WT transcripts",
     ylab="kb/min",
     lwd = 2,
     type = "l")

lines(rates_Maslon_m6A_WT$meRIP_WT,
      col="grey",
      lwd = 1,
      type = "l")
legend("topright", inset=.02, legend=c("Elongation rate", "m6A levels"),
       col=c("blue", "grey"), lty=1, lwd=2, cex=0.8)
dev.off()	

# ------------------
# Correlation test |
# ------------------
cor.test(saved_values_WT$Reads_norm_WT, saved_values_WT$elongation_b_min, method=c("pearson", "kendall", "spearman"))

# Pearson resulted to be more significant
png(file = paste0(rates_out, "5_Plots/Pearson_MaslonElongation_meRIP_WT_2Kb.png"))
ggscatter(saved_values_WT, y = "meRIP_WT", x = "Maslon_Kb_min",
          color="lightblue3",shape = 20,
          size=1.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "m6A levels", xlab = "Elongation rate Maslon data (Kb/min)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Pearson correlation m6A levels vs. elongation rates (Maslon data) - WT (746 genes + 2Kb)")
dev.off()
