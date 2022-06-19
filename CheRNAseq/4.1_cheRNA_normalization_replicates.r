#!/usr/bin/env Rscript

#########################################
## Normalize cheRNAs (RefSeq gene list) #
#########################################

# -------
# Paths |
# -------
cheRNA_path <- "/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/1_Intersect_RefSeq/"
output <- "/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/2_Normalized_data/"
output_plots <- "/media/cc/A/Alicia/NGS/cheRNA/cheRNA_output/2_Normalized_data/3_Plots/"

# ------------
# Open files |
# ------------
cheRNA_list = list.files(cheRNA_path, pattern="*.bed")
cheRNA_list

for (i in seq_along(cheRNA_list)) {
  filename <- sub("_RefSeq_intersect.bed", "", cheRNA_list[i])
  df <- read.table(paste0(cheRNA_path,cheRNA_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",
                                 "cheRNA"))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(filename, df)
}

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
TKO_I_reads	<-   93717203
TKO_II_reads	<- 93903878
TKO_III_reads	<- 98738497
WT_I_reads	<-   85262714
WT_II_reads	<-   98517689
WT_III_reads	<- 83632131

# ---------------
# Normalization |
# ---------------

# Normalization - total reads / transcript size 
WT_I$cheRNA_normI = ((WT_I$cheRNA/WT_I_reads) / ((WT_I$End-WT_I$Start)/1000)) 
WT_II$cheRNA_normII = ((WT_II$cheRNA/WT_II_reads) / ((WT_II$End-WT_II$Start)/1000)) 
WT_III$cheRNA_normIII = ((WT_III$cheRNA/WT_III_reads) / ((WT_III$End-WT_III$Start)/1000)) 

TKO_I$cheRNA_normI = ((TKO_I$cheRNA/TKO_I_reads) / ((TKO_I$End-TKO_I$Start)/1000)) 
TKO_II$cheRNA_normII = ((TKO_II$cheRNA/TKO_II_reads) / ((TKO_II$End-TKO_II$Start)/1000)) 
TKO_III$cheRNA_normIII = ((TKO_III$cheRNA/TKO_III_reads) / ((TKO_III$End-TKO_III$Start)/1000)) 

# ------------------------------
# Get means between replicates |
# ------------------------------
WT_table <- subset(WT_I, select = -c(cheRNA))
TKO_table <- subset(TKO_I, select = -c(cheRNA))

sum(WT_I$Gene_name==WT_III$Gene_name) # chack all rows are equal between tables
sum(TKO_I$Gene_name==TKO_III$Gene_name) # chack all rows are equal between tables

WT_table$cheRNA_normII   = WT_II$cheRNA_normII
WT_table$cheRNA_normIII   = WT_III$cheRNA_normIII
TKO_table$cheRNA_normII   = TKO_II$cheRNA_normII
TKO_table$cheRNA_normIII   = TKO_III$cheRNA_normIII

WT_table$cheRNA_norm_mean = (WT_table$cheRNA_normI + WT_table$cheRNA_normII + WT_table$cheRNA_normIII)/3
TKO_table$cheRNA_norm_mean = (TKO_table$cheRNA_normI + TKO_table$cheRNA_normII + TKO_table$cheRNA_normIII)/3

# ------------
# Merge data |
# ------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                  list(WT_table,
                       TKO_table[,c("cheRNA_normI","cheRNA_normII","cheRNA_normIII","cheRNA_norm_mean","Gene_name")]))

colnames(merged) <- c("Gene_name", "Chr","Start","End","Strand",
                            "cheRNA_WT1","cheRNA_WT2","cheRNA_WT3","cheRNA_WT_mean",
                            "cheRNA_TKO1","cheRNA_TKO2","cheRNA_TKO3","cheRNA_TKO_mean")

# -----------
# Save data |
# -----------
write.table(merged, 
            file = paste0(output,"cheRNA_RefSeq_LongList_Normalized.bed"),
            quote = F, sep="\t", col.names = T, row.names = F)     


# ----------
# Boxplots |
# ----------

png(file = paste0(output_plots, "BoxPlots_cheRNAs_replicates.png"))
boxplot(merged$cheRNA_WT1, merged$cheRNA_WT3, merged$cheRNA_WT3,
        merged$cheRNA_TKO1, merged$cheRNA_TKO2, merged$cheRNA_TKO3,
        col = c(4,4,4,2,2,2),
        at = c(1,2,3,5,6,7),
        names = c("WT1","WT2","WT3", "TKO1","TKO2","TKO3"),
        main = "cheRNA (replicates) (20621 genes)",
        ylab="cheRNA reads",
        outline = F,
        boxwex = 0.5)
dev.off()	

#wilcox <- wilcox.test(merged$cheRNA_WT_mean, merged$cheRNA_TKO_mean)
#wilcox
png(file = paste0(output_plots, "BoxPlots_cheRNAs_means.png"))
boxplot(merged$cheRNA_WT_mean, merged$cheRNA_TKO_mean,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "cheRNA (means) (20621 genes)",
        ylab="cheRNA reads",
        outline = F,
        xlab=paste("Wilcox test: W = 208182100, p-value = 0.0002448"),
        boxwex = 0.3)
dev.off()	

merged_ordered <- merged[order(merged$cheRNA_WT_mean, decreasing = T),]
merged_ordered2 <- merged[order(merged$cheRNA_TKO_mean, decreasing = T),]
plot(density(merged$cheRNA_WT1, to=0.000001), 
     col="blue", 
     lty=1, 
     main = paste("WT and TKO elongation rates ingenes (threshold=0.7) \n (genes with at least one NA were removed in all the samples)"), 
     cex.main=1, 
     xlab = "Elongation rate (kb/min) \n bandwidth = 0.6")
lines(density(merged$cheRNA_WT2, to=0.000001), lty=1, col="blue")
lines(density(merged$cheRNA_WT3, to=0.000001), lty=1, col="blue")
lines(density(merged$cheRNA_TKO1, to=0.000001), lty=1, col="red")
lines(density(merged$cheRNA_TKO2, to=0.000001), lty=1, col="red")
lines(density(merged$cheRNA_TKO3, to=0.000001), lty=1, col="red")

hist(merged$cheRNA_WT_mean,
     breaks=1000,
     xlim=c(0,0.000005),
     ylim=c(0,13000),
     col="blue",
     border='blue',
     #main = paste0(number," genes in WT"),
     xlab = paste("Elongation rate (kb/min) \n bin size = 500"))
hist(merged$cheRNA_TKO_mean,
     breaks=1000,
     xlim=c(0,0.000005),
     ylim=c(0,13000),
     col="red",
     border='red',
     #main = paste0(number," genes in WT"),
     xlab = paste("Elongation rate (kb/min) \n bin size = 500"))

# ----------------------------
# Mann-Whitney-Wilcoxon test |
# ----------------------------

# GB
sink(paste0(output_plots, "RNApolII_WilcoxTest_GB_500bp.txt"))
wilcox.test(mergedRNAPII_GB$RNAPII_WTmean_GB, mergedRNAPII_GB$RNAPII_TKOmean_GB)
sink()

# PR
sink(paste0(output_plots, "RNApolII_WilcoxTest_PR_500bp.txt"))
wilcox.test(mergedRNAPII_PR$RNAPII_WTmean_PR, mergedRNAPII_PR$RNAPII_TKOmean_PR)
sink()