#!/usr/bin/env Rscript

##########################################
## Normalize RNApolII (RefSeq long list) #
##########################################

# -------
# Paths |
# -------
RNAPII_GB_path <- "/path/RNAPII/Intersect_bedtools/GeneBody/"
RNAPII_PR_path <- "/path/RNAPII/Intersect_bedtools/Promoters/"

output <- "/path/RNAPII/Normalized_data/"
output_plots <- "/path/RNAPII/BoxPlots/"

# -------------------------------------
# Open BED data meRIP-Ensembl aligned |
# -------------------------------------
GB_list = list.files(RNAPII_GB_path, pattern="*.bed")
PR_list = list.files(RNAPII_PR_path, pattern="*.bed")

for (i in seq_along(GB_list)) {
  filename <- sub("_.*", "", GB_list[i])
  df <- read.table(paste0(RNAPII_GB_path,GB_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",
                                 "RNAPII_GB"))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(paste0(filename,"_GB"), df)
}

for (i in seq_along(PR_list)) {
  filename <- sub("_.*", "", PR_list[i])
  df <- read.table(paste0(RNAPII_PR_path,PR_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",
                                 "RNAPII_PR"))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(paste0(filename,"_PR"), df)
}

# ------------------------------------
# Total reads (data from experiment) |
# ------------------------------------
WU1_reads	<- 58416135
WU2_reads	<- 62366120
TU1_reads	<- 50034143
TU2_reads	<- 51396457

# ----------------------
# Nomalization factors |
# ----------------------
WT_factor  <- 1
TKO_factor <- 1.155026393

# ---------------
# Normalization |
# ---------------

# GB normalization - total reads / transcript size / normalization factor
WT1_GB$RNAPII_GB_norm = ((WT1_GB$RNAPII_GB/WU1_reads) / ((WT1_GB$End-WT1_GB$Start)/1000)) * WT_factor 
WT2_GB$RNAPII_GB_norm = ((WT2_GB$RNAPII_GB/WU2_reads) / ((WT2_GB$End-WT2_GB$Start)/1000)) * WT_factor 
TKO1_GB$RNAPII_GB_norm = ((TKO1_GB$RNAPII_GB/TU1_reads) / ((TKO1_GB$End-TKO1_GB$Start)/1000)) * TKO_factor
TKO2_GB$RNAPII_GB_norm = ((TKO2_GB$RNAPII_GB/TU2_reads) / ((TKO2_GB$End-TKO2_GB$Start)/1000)) * TKO_factor 

# PR normalization - total reads / normalization factor
WT1_PR$RNAPII_PR_norm = (WT1_PR$RNAPII_PR/WU1_reads) * WT_factor 
WT2_PR$RNAPII_PR_norm = (WT2_PR$RNAPII_PR/WU2_reads) * WT_factor 
TKO1_PR$RNAPII_PR_norm = (TKO1_PR$RNAPII_PR/TU1_reads) * TKO_factor 
TKO2_PR$RNAPII_PR_norm = (TKO2_PR$RNAPII_PR/TU2_reads) * TKO_factor 

# ------------------------------
# Get means between replicates |
# ------------------------------
WT_GB <- subset(WT1_GB, select = -c(RNAPII_GB))
TKO_GB <- subset(TKO1_GB, select = -c(RNAPII_GB))
WT_PR <- subset(WT1_PR, select = -c(RNAPII_PR))
TKO_PR <- subset(TKO1_PR, select = -c(RNAPII_PR))

sum(WT2_GB$Gene_name==TKO2_GB$Gene_name) # chack all rows are equal between tables
sum(WT1_PR$Gene_name==WT2_PR$Gene_name) # chack all rows are equal between tables

WT_GB$RNAPII_GB_WT2   = WT2_GB$RNAPII_GB_norm
TKO_GB$RNAPII_GB_TKO2 = TKO2_GB$RNAPII_GB_norm
WT_PR$RNAPII_PR_WT2   = WT2_PR$RNAPII_PR_norm
TKO_PR$RNAPII_PR_TKO2 = TKO2_PR$RNAPII_PR_norm

WT_GB$RNAPII_GB_WT = (WT_GB$RNAPII_GB_norm+WT_GB$RNAPII_GB_WT2)/2
TKO_GB$RNAPII_GB_TKO = (TKO_GB$RNAPII_GB_norm+TKO_GB$RNAPII_GB_TKO2)/2
WT_PR$RNAPII_PR_WT = (WT_PR$RNAPII_PR_norm+WT_PR$RNAPII_PR_WT2)/2
TKO_PR$RNAPII_PR_TKO = (TKO_PR$RNAPII_PR_norm+TKO_PR$RNAPII_PR_TKO2)/2

# ------------
# Merge data |
# ------------
mergedRNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                  list(WT_GB,
                       TKO_GB[,c("RNAPII_GB_norm","RNAPII_GB_TKO2","RNAPII_GB_TKO","Gene_name")]))

colnames(mergedRNAPII_GB) <- c("Gene_name", "Chr","Start","End","Strand",
                            "RNAPII_WT1_GB","RNAPII_WT2_GB","RNAPII_WTmean_GB",
                            "RNAPII_TKO1_GB","RNAPII_TKO2_GB","RNAPII_TKOmean_GB")

mergedRNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                         list(WT_PR,
                              TKO_PR[,c("RNAPII_PR_norm", "RNAPII_PR_TKO2", "RNAPII_PR_TKO", "Gene_name")]))

colnames(mergedRNAPII_PR) <- c("Gene_name","Chr","Start","End", "Strand",
                            "RNAPII_WT1_PR","RNAPII_WT2_PR","RNAPII_WTmean_PR",
                            "RNAPII_TKO1_PR","RNAPII_TKO2_PR","RNAPII_TKOmean_PR")

# -----------
# Save data |
# -----------
write.table(mergedRNAPII_GB, 
            file = paste0(output,"GeneBody/RNAPII_RefSeq_LongList_GB_500bp.bed"),
            quote = F, sep="\t", col.names = T, row.names = F)     

write.table(mergedRNAPII_PR, 
            file = paste0(output,"Promoters/RNAPII_RefSeq_LongList_PR_500bp.bed"),
            quote = F, sep="\t", col.names = T, row.names = F)     

# ----------
# Boxplots |
# ----------

# GB
png(file = paste0(output_plots, "BoxPlots_RNApolII_RefSeqLongList_GB_replicates_500bp.png"))
boxplot(mergedRNAPII_GB$RNAPII_WT1_GB, mergedRNAPII_GB$RNAPII_WT2_GB,
        mergedRNAPII_GB$RNAPII_TKO1_GB, mergedRNAPII_GB$RNAPII_TKO2_GB,
        col = c(4,4,2,2),
        at = c(1,2,4,5),
        names = c("WT1","WT2", "TKO1","TKO2"),
        main = "RNApolII Gene Body (replicates) (20546 genes)",
        ylab="RNApolII occupancy",
        outline = F,
        boxwex = 0.5)
dev.off()	

png(file = paste0(output_plots, "BoxPlots_RNApolII_RefSeqLongList_GB_means_500bp.png"))
boxplot(mergedRNAPII_GB$RNAPII_WTmean_GB, mergedRNAPII_GB$RNAPII_TKOmean_GB,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "RNApolII Gene Body (means) (20546 genes)",
        ylab="RNApolII occupancy",
        outline = F,
        xlab=paste("Wilcox test: W = 174731727, p-value < 2.2e-16"),
        boxwex = 0.3)
dev.off()	

# PR
png(file = paste0(output_plots, "BoxPlots_RNApolII_RefSeqLongList_PR_replicates_500bp.png"))
boxplot(mergedRNAPII_PR$RNAPII_WT1_PR, mergedRNAPII_PR$RNAPII_WT2_PR,
        mergedRNAPII_PR$RNAPII_TKO1_PR, mergedRNAPII_PR$RNAPII_TKO2_PR,
        col = c(4,4,2,2),
        at = c(1,2,4,5),
        names = c("WT1","WT2", "TKO1","TKO2"),
        main = "RNApolII Promoters (replicates) (20546 genes)",
        ylab="RNApolII occupancy",
        outline = F,
        boxwex = 0.5)
dev.off()	

png(file = paste0(output_plots, "BoxPlots_RNApolII_RefSeqLongList_PR_means_500bp.png"))
boxplot(mergedRNAPII_PR$RNAPII_WTmean_PR, mergedRNAPII_PR$RNAPII_TKOmean_PR,
        col = c(4,2),
        at = c(1,2),
        names = c("WT", "TKO"),
        main = "RNApolII Promoters (means) (20546 genes)",
        ylab="RNApolII occupancy",
        outline = F,
        xlab=paste("Wilcox test: W = 194923900, p-value < 2.2e-16"),
        boxwex = 0.3)
dev.off()	

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
