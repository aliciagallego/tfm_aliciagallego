#!/usr/bin/env Rscript

#############################################
## Pausing index (RefSeq long list -+500bp) #
#############################################

# -----------
# Libraries |
# -----------
library("ggpubr")

# -------
# Paths |
# -------
RNAPII_GB_path <- "/path/RNAPII/Intersect_bedtools/GeneBody/"
RNAPII_PR_path <- "/path/RNAPII/Intersect_bedtools/Promoters/"

m6AWT_path <- "/path/meRIP/Normalized_data/meRIP_RefSeqLongList_Normalized_WT_2Kb.txt"
m6ATKO_path <-"/path/meRIP/Normalized_data/meRIP_RefSeqLongList_Normalized_TKO_2Kb.txt"

output <- "/path/RNAPII/Pausing_index/"
output_plots <- "/path/RNAPII/Pausing_index/Plots/"
output_correlation <- "/path/RNAPII/Pausing_index/Correlations/PausingIndex_m6A/"

# -----------
# Open data |
# -----------
RNAPII_GB <- read.table(RNAPII_GB_path ,h=T, sep="\t", stringsAsFactors=F)
RNAPII_PR <- read.table(RNAPII_PR_path ,h=T, sep="\t", stringsAsFactors=F)

m6AWT <- read.table(m6AWT_path, h=T,sep="\t",stringsAsFactors=F)
m6ATKO <- read.table(m6ATKO_path, h=T,sep="\t",stringsAsFactors=F)

# ---------------------
# Merge RNApolII data |
# ---------------------
mergedRNAPII_GB_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                         list(RNAPII_GB[,c("Gene_name","Chr","Start","End","Strand","RNAPII_WTmean_GB","RNAPII_TKOmean_GB")],
                              RNAPII_PR[,c("Gene_name","Start","End","RNAPII_WTmean_PR","RNAPII_TKOmean_PR")]))

colnames(mergedRNAPII_GB_PR) <- c("Gene_name","Chr","Start_GB","End_GB","Strand","RNAPII_WT_GB","RNAPII_TKO_GB",
                               "Start_PR","End_PR","RNAPII_WT_PR","RNAPII_TKO_PR")

# -------------------------
# Calculate pausing index |
# -------------------------
mergedRNAPII_GB_PR$PausIndex_WT = mergedRNAPII_GB_PR$RNAPII_WT_PR/mergedRNAPII_GB_PR$RNAPII_WT_GB
mergedRNAPII_GB_PR$PausIndex_TKO = mergedRNAPII_GB_PR$RNAPII_TKO_PR/mergedRNAPII_GB_PR$RNAPII_TKO_GB

# -----------------------------
# Remove NA, Inf and 0 values |
# -----------------------------
mergedRNAPII_GB_PR <- subset(mergedRNAPII_GB_PR, PausIndex_WT!="NaN" & PausIndex_WT!=0 & PausIndex_WT!="Inf")
mergedRNAPII_GB_PR <- subset(mergedRNAPII_GB_PR, PausIndex_TKO!="NaN" & PausIndex_TKO!=0 & PausIndex_TKO!="Inf")

# -----------
# Save data |
# -----------
write.table(mergedRNAPII_GB_PR, 
            file = paste0(output,"RNAPII_RefSeq_LongList_PausingIndex_500bp.txt"),
            quote = F, sep="\t", col.names = T, row.names = F)     

# ----------------------------------------
# Merge RNAPII pausing data and m6A data |
# ----------------------------------------
merged<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                     m6ATKO[,c("meRIP_TKO","Gene_name")],
                     mergedRNAPII_GB_PR[,c("PausIndex_WT","PausIndex_TKO", "Gene_name")]))

# ---------------------------------------------------------------
# Remove outliers from meRIP data and log transform RNAPII data |
# ---------------------------------------------------------------
Q <- quantile(merged$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(merged, merged$meRIP_WT > low & merged$meRIP_WT < up)
saved_values_WT$PausIndex_WT <- log(saved_values_WT$PausIndex_WT+1)

Q <- quantile(merged$meRIP_TKO, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(merged$meRIP_TKO)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_TKO <- subset(merged, merged$meRIP_TKO > low & merged$meRIP_TKO < up)
saved_values_TKO$PausIndex_TKO <- log(saved_values_TKO$PausIndex_TKO+1)

# ---------------------------------
# Boxplots RNApolII pausing index |
# ---------------------------------
meanwT <- round(mean(mergedRNAPII_GB_PR$PausIndex_WT),digits=2)
meanTKO <- round(mean(mergedRNAPII_GB_PR$PausIndex_TKO), digits=2)
means <- c(meanwT,meanTKO)

png(file = paste0(output_plots, "BoxPlots_RNAPII_RefSeq_LongList_PausingIndex_500bp.png"))
boxplot(mergedRNAPII_GB_PR$PausIndex_WT, mergedRNAPII_GB_PR$PausIndex_TKO,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "RNApolII pausing index (means) (20217 genes)",
        ylab="Promoter pausing index",
        outline = F,
        xlab=paste("Wilcox test: W = 207061353, p-value = 0.02151"),
        boxwex = 0.3)
text(x=c(1:2),
     y=4.2,
     paste("mean = ",means,sep=""))
dev.off()	

# ---------------------------------------------------
# Mann-Whitney-Wilcoxon test RNApolII pausing index |
# ---------------------------------------------------
sink(paste0(output_plots, "RNApolII_WilcoxTest_RefSeq_LongList_PausingIndex_500bp.txt"))
wilcox.test(mergedRNAPII_GB_PR$PausIndex_WT, mergedRNAPII_GB_PR$PausIndex_TKO)
sink()

# -----------------------------------------------
# Correlation test RNAPII pausing index and m6A |
# -----------------------------------------------
genesWT <- nrow(saved_values_WT)
genesTKO <- nrow(saved_values_TKO)

png(file = paste0(output_correlation, "Pearson_PausingIndex_meRIP_WT_500bp_2.png"))
ggscatter(saved_values_WT, x = "meRIP_WT", y = "PausIndex_WT",
          color="steelblue1",shape = 4,
          size=0.5,
          add.params = list(color = "steelblue3", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylim=c(-2,5),
          xlab = "m6A levels", ylab = "log RNApolII pausing index") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Pearson correlation RNApolII pausing index vs. m6A levels \n WT (",genesWT," genes, RefSeq Long List)"))
dev.off()

png(file = paste0(output_correlation, "Pearson_PausingIndex_meRIP_TKO_500bp_2.png"))
ggscatter(saved_values_TKO, x = "meRIP_TKO", y = "PausIndex_TKO",
          color="red",shape = 4,
          size=0.5,
          add.params = list(color = "red3", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylim=c(-2,5),
          xlab = "m6A levels", ylab = "log RNApolII pausing index") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Pearson correlation RNApolII pausing index vs. m6A levels \n TKO (",genesTKO," genes, RefSeq Long List)"))
dev.off()
