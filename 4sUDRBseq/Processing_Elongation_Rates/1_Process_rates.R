#!/usr/bin/env Rscript

###########################################
## Normalize elongation rate calculations # 
###########################################

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

output <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/")

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
                                TTseq[,c("Gene_name","WT_Pull", "TK0_Pull")]))

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

# ------------------
# Save merged file |
# ------------------
write.table(merged_non05, file = paste0(output,"1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull_processed_without05.txt"), 
            quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

# ----------------------------
# Merge rates and other data |
# ----------------------------
merged_WT_alldata<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                           list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                                RNApolGB[,c("RNAPII_WTmean_GB", "Gene_name")],
                                RNApolPR[,c("RNAPII_WTmean_PR", "Gene_name")],
                                Pausing[,c("PausIndex_WT","Gene_name")],
                                cheRNA[,c("cheRNA_WT_mean","Gene_name")],
                                merged_non05[,c("WT_Pull", "Exon_number", "Size","Gene_name")]))

merged_TKO_alldata<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                           list(m6ATKO[,c("Chr","Start","End","Strand","meRIP_TKO","Gene_name")],
                                RNApolGB[,c("RNAPII_TKOmean_GB", "Gene_name")],
                                RNApolPR[,c("RNAPII_TKOmean_PR", "Gene_name")],
                                Pausing[,c("PausIndex_TKO","Gene_name")],
                                cheRNA[,c("cheRNA_TKO_mean","Gene_name")],
                                merged_non05[,c("TK0_Pull","Exon_number", "Size","Gene_name")]))

# ----------------------------------------------
# Bind WT and TKO rows in one single dataframe |
# ----------------------------------------------
names_vec <- c("Gene_name", "Chr", "Start", "End", "Strand", 
               "meRIP","RNAPII_GB", "RNAPII_PR", "PausIndex", 
               "cheRNA", "Elongation_rate","Exon_number", "Size")

colnames(merged_WT_alldata) <- names_vec
colnames(merged_TKO_alldata) <- names_vec

merged_WT_alldata$label <- "0"
merged_TKO_alldata$label <- "1"

length(merged_TKO_alldata$label)
bind_WT_TKO <- rbind(merged_WT_alldata, merged_TKO_alldata)

# -------------------
# Save merged files |
# -------------------
write.table(merged_WT_alldata, file = paste0(output,"7_Rates_and_all_data/Rates_and_all_data_without05_WT.txt"), 
            quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

write.table(merged_TKO_alldata, file = paste0(output,"7_Rates_and_all_data/Rates_and_all_data_without05_TKO.txt"), 
            quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

write.table(bind_WT_TKO, file = paste0(output,"7_Rates_and_all_data/Rates_and_all_data_without05_WT_TKO_bind.txt"), 
            quote = F, sep="\t", col.names = TRUE, row.names = FALSE)

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
nrow(TKOfast)

# m6A
merged_WTslow_m6A<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                                WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTmedium_m6A<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                                WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTfast_m6A<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(m6AWT[,c("Chr","Start","End","Strand","meRIP_WT","Gene_name")],
                                WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_TKOslow_m6A<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                           list(m6ATKO[,c("Chr","Start","End","Strand","meRIP_TKO","Gene_name")],
                                TKOslow[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOmedium_m6A<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                            list(m6ATKO[,c("Chr","Start","End","Strand","meRIP_TKO","Gene_name")],
                                 TKOmedium[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOfast_m6A<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                            list(m6ATKO[,c("Chr","Start","End","Strand","meRIP_TKO","Gene_name")],
                                 TKOfast[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# RNAPII WT
merged_WTslow_RNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_WTmean_GB","Gene_name")],
                                      WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTmedium_RNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_WTmean_GB", "Gene_name")],
                                      WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTfast_RNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                   list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_WTmean_GB", "Gene_name")],
                                        WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTslow_RNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                   list(RNApolPR[,c("Chr","Start","End","Strand","RNAPII_WTmean_PR","Gene_name")],
                                        WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTmedium_RNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(RNApolPR[,c("Chr","Start","End","Strand","RNAPII_WTmean_PR","Gene_name")],
                                      WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTfast_RNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(RNApolPR[,c("Chr","Start","End","Strand","RNAPII_WTmean_PR","Gene_name")],
                                      WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

# RNAPII TKO
merged_TKOslow_RNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                 list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_TKOmean_GB", "Gene_name")],
                                      TKOslow[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOmedium_RNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                  list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_TKOmean_GB", "Gene_name")],
                                       TKOmedium[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOfast_RNAPII_GB<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                  list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_TKOmean_GB", "Gene_name")],
                                       TKOfast[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOslow_RNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                  list(RNApolPR[,c("Chr","Start","End","Strand","RNAPII_TKOmean_PR","Gene_name")],
                                       TKOslow[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOmedium_RNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                  list(RNApolPR[,c("Chr","Start","End","Strand","RNAPII_TKOmean_PR","Gene_name")],
                                       TKOmedium[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKO_fast_RNAPII_PR<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                  list(RNApolPR[,c("Chr","Start","End","Strand","RNAPII_TKOmean_PR","Gene_name")],
                                       TKOfast[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# Pausing index
merged_WTslow_Pausing<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                               list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_WT","Gene_name")],
                                    WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTmedium_Pausing<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                               list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_WT","Gene_name")],
                                    WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTfast_Pausing<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                               list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_WT","Gene_name")],
                                    WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_TKOslow_Pausing<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                               list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_TKO","Gene_name")],
                                    TKOslow[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOmedium_Pausing<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_TKO","Gene_name")],
                                     TKOmedium[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOfast_Pausing<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                                list(Pausing[,c("Chr","Start_GB","End_GB","Strand","PausIndex_TKO","Gene_name")],
                                     TKOfast[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# cheRNAs
merged_WTslow_cheRNA<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                              list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_WT_mean","Gene_name")],
                                   WTslow[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTmedium_cheRNA<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                              list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_WT_mean","Gene_name")],
                                   WTmedium[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_WTfast_cheRNA<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                              list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_WT_mean","Gene_name")],
                                   WTfast[,c("TTseq_WT_pull", "Exon_number","Size","Gene_name")]))

merged_TKOslow_cheRNA<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                              list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_TKO_mean","Gene_name")],
                                   TKOslow[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOmedium_cheRNA<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                               list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_TKO_mean","Gene_name")],
                                    TKOmedium[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

merged_TKOfast_cheRNA<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                               list(cheRNA[,c("Chr","Start","End","Strand","cheRNA_TKO_mean","Gene_name")],
                                    TKOfast[,c("TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# ---------------------------
# Identifying biased values |
# ---------------------------
number <- nrow(merged_non05)

png(file = paste0(output, "5_Plots_without05/Histograms_Pull_WT_rates.png"))
hist(merged_non05$WT_Pull,
     breaks=500,
     col="blue",
     border='blue',
     main = paste0(number," genes in WT"),
     xlab = paste("Elongation rate (kb/min) \n bin size = 500"))
dev.off()

png(file = paste0(output, "5_Plots_without05/Histograms_Pull_TKO_rates.png"))
hist(merged_non05$TK0_Pull,
     breaks=500,
     col="red",
     border='red',
     main = paste0(number," genes in TKO"),
     xlab = paste("Elongation rate (kb/min) \n bin size = 500"))
dev.off()

# -------------------------------
# Elongation rate Density plots |
# -------------------------------
number_non05 <- nrow(merged_non05)

png(file = paste0(output, "5_Plots_without05/Density_Pull_NAall_bw02_higher05values.png"))
plot(density(merged_non05$WT_Pull,to=7, bw=0.2), 
     col="blue",  
     main = paste("WT and TKO elongation rates in", number_non05, "genes"), 
     cex.main=1,
     lwd=2,
     xlab = "Elongation rate (kb/min) \n bandwidth = 0.2")
lines(density(merged_non05$TK0_Pull, to=7, bw=0.2), col="red",lwd=2)
legend("topright", c("WT","TKO"), 
       col=c("blue","red"), lwd=2, inset = .02,cex=0.9)
dev.off()

png(file = paste0(output, "5_Plots_without05/Density_Pull_NAall_bw06_higher05values.png"))
plot(density(merged_non05$WT_Pull,to=7, bw=0.6), 
     col="blue",  
     main = paste("WT and TKO elongation rates in", number_non05, "genes (threshold=0.7) \n (only genes with rates > 0.5 Kb/min)"), 
     cex.main=1, 
     xlab = "Elongation rate (kb/min) \n bandwidth = 0.6")
lines(density(merged_non05$TK0_Pull, to=7, bw=0.6), col="red")
legend("topright", c("WT","TKO"), 
       col=c("blue","red"), lwd=1.5, inset = .02,cex=0.9)
dev.off()

# -----------------------
# Density per quantiles |
# -----------------------

# bandwidth = 0.2
png(file = paste0(output, "5_Plots_without05/Density_Pull_tertiles_bw02.png"))
par(cex.lab=1) # is for y-axis
par(cex.axis=1.2) # is for x-axis
plot(density(WTmedium$TTseq_WT_pull, from=0, to=6, bw=.2), 
     mgp=c(3,0.5,0),
     col="blue2", 
     lty=2, 
     lwd=2,
     cex.lab=1.2,
     main = paste("WT and TKO elongation rates in", number_non05, "genes"), 
     cex.main=1, 
     xlab = "Elongation rate (kb/min) \n bandwidth = 0.2")
lines(density(WTslow$TTseq_WT_pull, from=0, to=8, bw=.2), lty=3, lwd=2, col="blue")
lines(density(WTfast$TTseq_WT_pull, from=0, to=8, bw=.2), lwd=2, col="blue3", lty=1)

lines(density(TKOslow$TTseq_TKO_pull, from=0, to=8, bw=.2), lwd=2, col="red", lty=3)
lines(density(TKOmedium$TTseq_TKO_pull, from=0, to=8, bw=.2), lwd=2, lty=2, col="red2")
lines(density(TKOfast$TTseq_TKO_pull, from=0, to=8, bw=.2), lwd=2, lty=1, col="red3")

legend("topright", c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"), 
       col=c("blue","blue2","blue3","red","red2","red3"), lwd=2, lty=c(3,2,1), inset = .02,cex=0.9)
dev.off()	

# bandwidth = 0.6
png(file = paste0(output, "5_Plots_without05/Density_Pull_tertiles_bw06.png"))
plot(density(WTmedium$WT_Pull, from=-1, to=8, bw=.6), 
     col="deepskyblue3", 
     lty=2, 
     main = paste("WT and TKO elongation rates in", number, "genes (threshold=0.7) \n (genes with at least one NA were removed in all the samples)"), 
     cex.main=1, 
     xlab = "Elongation rate (kb/min) \n bandwidth = 0.6")
lines(density(WTslow$WT_Pull, to=8, bw=.6), lty=3, col="deepskyblue")
lines(density(WTfast$WT_Pull, to=8, bw=.6), col="deepskyblue4", lty=1)

lines(density(TKOslow$TK0_Pull, to=8, bw=.6), col="indianred1", lty=3)
lines(density(TKOmedium$TK0_Pull, to=8, bw=.6), lty=2, col="indianred4")
lines(density(TKOfast$TK0_Pull, to=8, bw=.6), lty=1, col="red3")

legend("topright", c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"), 
       col=c("deepskyblue","deepskyblue3","deepskyblue4","indianred1","indianred4","red3"), lwd=1.5, lty=c(3,2,1), inset = .02,cex=0.9)
dev.off()	

# --------------------------
# Elongation rate Boxplots |
# --------------------------
row_noNA <- nrow(merged_non05)

png(file = paste0(output, "5_Plots_without05/Boxplots_Pull.png"))
wilcox <- wilcox.test(merged_non05$WT_Pull, merged_non05$TK0_Pull)
wilcox
par(cex.lab=1.5) # is for y-axis
par(cex.axis=1.5) # is for x-axis
boxplot(merged_non05$WT_Pull, merged_non05$TK0_Pull,
        col = c(4,2),
        at = c(1,2),
        names = c("WT","H1-TKO"),
        ylab="Elongation rates (kb/min)",
        main = paste0("WT and H1-TKO elongation rates in ",row_noNA," genes"),
        ylim=c(0,7),
        outline = F,
        cex.xlab=1.8,
        boxwex = 0.3)
#legend("topleft", c("WT","H1-TKO"), col=c("blue","red"), lwd=4, cex=1.5, inset = .02)
dev.off()	

png(file = paste0(output, "5_Plots_without05/Boxplots_Pull_tertiles_signif.png"))
#wilcoxfast <- wilcox.test(WTfast$WT_Pull, TKOfast$TK0_Pull)
#wilcoxfast
par(cex.lab=1.2) # is for y-axis
par(cex.axis=1.2) # is for x-axis
par(mar=c(5,4,4,2)+0.5)
boxplot(WTslow$TTseq_WT_pull,  TKOslow$TTseq_TKO_pull,  
        WTmedium$TTseq_WT_pull, TKOmedium$TTseq_TKO_pull, 
        WTfast$TTseq_WT_pull, TKOfast$TTseq_TKO_pull,
        col = c("blue","red","blue","red","blue","red"),
        at = c(1,2,3,4,5,6),
        names = c("WT slow","TKO slow","WT med","TKO med","WT fast","TKO fast"),
        ylab="Elongation rate (kb/min)",
        main = "WT and TKO elongation rates - means",
        las=2,
        ylim=c(0,6),
        cex.lab=1,
        outline = F,
        boxwex = 0.6)
segments(3,3.8,4,3.8,lwd=0.8)
segments(5,5.5,6,5.5,lwd=0.8)
text(3.5, 4, "*",cex=2,col=("black"))
text(5.5, 5.8, "*",cex=2,col=("black"))
legend("topleft", c("WT","H1-TKO"), col=c("blue","red"), lwd=4, inset = .02)
legend("bottomright", c("slow p = 9.22e-10","med p* < 2.2e-16", "fast  p* < 2.2e-16"), inset = .02, cex=1.5)
dev.off()	

# Tertiles & m6A
png(file = paste0(output, "5_Plots_without05/Boxplots_Pull_tertiles_m6A_signif.png"))
#wilcoxfast <- wilcox.test(merged_WTfast_m6A$meRIP_WT,  merged_TKOfast_m6A$meRIP_TKO)
#wilcoxfast
boxplot(merged_WTslow_m6A$meRIP_WT, merged_WTmedium_m6A$meRIP_WT, merged_WTfast_m6A$meRIP_WT,  
        merged_TKOslow_m6A$meRIP_TKO, merged_TKOmedium_m6A$meRIP_TKO, merged_TKOfast_m6A$meRIP_TKO,
        col = c("blue","blue2","blue3","red","red2","red3"),
        at = c(1,2,3,5,6,7),
        names = c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"),
        ylab="m6A levels",
        main = "m6A levels across elongation rate groups",
        ylim=c(-0.2,1.5),
        las=2,
        outline = F,
        boxwex = 0.6)
segments(1,1.26,5,1.26,lwd=0.8)
segments(2,1.4,6,1.4,lwd=0.8)
text(3, 1.3,"*",cex=2,col=("black"))
text(4, 1.44, "*",cex=2,col=("black"))
legend("bottomright", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
dev.off()	

# Tertiles & RNAPII_GB
par(mar=c(5,5,4,2)+0.6)
png(file = paste0(output, "5_Plots_without05/Boxplots_Pull_tertiles_RNAPII_GB.png"))
#wilcoxfast <- wilcox.test(merged_TKOmedium_RNAPII_GB$RNAPII_TKOmean_GB, merged_TKOfast_RNAPII_GB$RNAPII_TKOmean_GB)
#wilcoxfast
boxplot(merged_WTslow_RNAPII_GB$RNAPII_WTmean_GB, merged_WTmedium_RNAPII_GB$RNAPII_WTmean_GB, merged_WTfast_RNAPII_GB$RNAPII_WTmean_GB,
        merged_TKOslow_RNAPII_GB$RNAPII_TKOmean_GB, merged_TKOmedium_RNAPII_GB$RNAPII_TKOmean_GB, merged_TKOfast_RNAPII_GB$RNAPII_TKOmean_GB,
        col = c("blue","blue2","blue3","red","red2","red3"),
        at = c(1,2,3,5,6,7),
        names = c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"),
        #ylab="RNAPII occupancies in Gene Body",
        main = "RNAPII in gene bodies\nacross elongation rate groups",
        ylim=c(0,1.7e-06),
        las=2,

        outline = F,
        boxwex = 0.6)
segments(1,1.2e-06,2,1.2e-06,lwd=0.8)
segments(1,1.4e-06,3,1.4e-06,lwd=0.8)
segments(5,1.3e-06,6,1.3e-06,lwd=0.8)
segments(5,1.45e-06,7,1.45e-06,lwd=0.8)

text(1.5, 1.25e-06,"*",cex=2,col=("black"))
text(2, 1.45e-06, "*",cex=2,col=("black"))
text(5.5, 1.35e-06,"*",cex=2,col=("black"))
text(6, 1.5e-06, "*",cex=2,col=("black"))

legend("topleft", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
legend("bottomright", c("slow-med p* < 2.2e-16","slow-fast  p* < 2.2e-16"), inset = .02, cex=.7)
dev.off()	

# Tertiles & RNAPII_PR
png(file = paste0(output, "5_Plots_without05/Boxplots_Pull_tertiles_RNAPII_PR.png"))
wilcoxfast <- wilcox.test(merged_WTfast_RNAPII_PR$RNAPII_WTmean_PR, merged_TKO_fast_RNAPII_PR$RNAPII_TKOmean_PR)
wilcoxfast
boxplot(merged_WTslow_RNAPII_PR$RNAPII_WTmean_PR, merged_WTmedium_RNAPII_PR$RNAPII_WTmean_PR, merged_WTfast_RNAPII_PR$RNAPII_WTmean_PR,
        merged_TKOslow_RNAPII_PR$RNAPII_TKOmean_PR, merged_TKOmedium_RNAPII_PR$RNAPII_TKOmean_PR, merged_TKO_fast_RNAPII_PR$RNAPII_TKOmean_PR,
        col = c("blue","blue2","blue3","red","red2","red3"),
        at = c(1,2,3,5,6,7),
        names = c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"),
        #ylab="RNAPII occupancies in promoters\nacross elongation rate groups",
        main = "RNAPII in promoters\nacross elongation rate groups",
        ylim=c(-0.6e-06,9e-06),
        las=2,
        outline = F,
        boxwex = 0.6)

segments(1,7e-06,2,7e-06,lwd=0.8)
segments(1,8e-06,3,8e-06,lwd=0.8)
segments(5,7.3e-06,6,7.3e-06,lwd=0.8)
segments(5,8e-06,7,8e-06,lwd=0.8)

text(1.5, 7.26e-06,"*",cex=2,col=("black"))
text(2, 8.26e-06, "*",cex=2,col=("black"))
text(5.5, 7.56e-06,"*",cex=2,col=("black"))
text(6, 8.26e-06, "*",cex=2,col=("black"))

legend("topleft", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
legend("bottomright", c("slow-med p* < 2.2e-16","slow-fast  p* < 2.2e-16"), inset = .02, cex=.7)
dev.off()	

# Tertiles & Pausing index
png(file = paste0(output, "5_Plots_without05/Boxplots_Pull_tertiles_PausingIndex.png"))
wilcoxfast <- wilcox.test(merged_WTfast_Pausing$PausIndex_WT,  merged_TKOfast_Pausing$PausIndex_TKO)
wilcoxfast
boxplot(merged_WTslow_Pausing$PausIndex_WT, merged_WTmedium_Pausing$PausIndex_WT, merged_WTfast_Pausing$PausIndex_WT,
        merged_TKOslow_Pausing$PausIndex_TKO, merged_TKOmedium_Pausing$PausIndex_TKO, merged_TKOfast_Pausing$PausIndex_TKO,
        col = c("blue","blue2","blue3","red","red2","red3"),
        at = c(1,2,3,5,6,7),
        names = c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"),
        ylab="Pausing index",
        main = "Pausing index\nacross elongation rate groups",
        ylim=c(-2,12),
        las=2,
        outline = F,
        boxwex = 0.6)
segments(1,10,2,10,lwd=0.8)
segments(1,11,3,11,lwd=0.8)
segments(5,10,6,10,lwd=0.8)
segments(5,11,7,11,lwd=0.8)

text(1.5, 10.5,"*",cex=2,col=("black"))
text(2, 11.5, "*",cex=2,col=("black"))
text(5.5, 10.5,"*",cex=2,col=("black"))
text(6, 11.5, "*",cex=2,col=("black"))

legend("topleft", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
legend("bottomright", c("slow-med p* < 2.2e-16","slow-fast  p* < 2.2e-16"), inset = .02, cex=.7)
dev.off()	

# Tertiles & cheRNAs
png(file = paste0(output, "5_Plots_without05/Boxplots_Pull_tertiles_cheRNAs.png"))
wilcoxfast <- wilcox.test(merged_WTfast_cheRNA$cheRNA_WT_mean,  merged_TKOmedium_cheRNA$cheRNA_TKO_mean)
wilcoxfast
boxplot(merged_WTslow_cheRNA$cheRNA_WT_mean, merged_WTmedium_cheRNA$cheRNA_WT_mean, merged_WTfast_cheRNA$cheRNA_WT_mean,
        merged_TKOslow_cheRNA$cheRNA_TKO_mean, merged_TKOmedium_cheRNA$cheRNA_TKO_mean, merged_TKOfast_cheRNA$cheRNA_TKO_mean,
        col = c("blue","blue2","blue3","red","red2","red3"),
        at = c(1,2,3,5,6,7),
        names = c("WT slow","WT med","WT fast","TKO slow","TKO med","TKO fast"),
        #ylab="cheRNA #reads",
        main = "cheRNA expression levels\nacross elongation rate groups",
        ylim=c(-1e-6,6e-06),
        las=2,
        outline = F,
        boxwex = 0.6)
segments(1,3.6e-06,2,3.6e-06,lwd=0.8)
segments(1,4.6e-06,3,4.6e-06,lwd=0.8)
segments(5,3.3e-06,6,3.3e-06,lwd=0.8)
segments(5,4.3e-06,7,4.3e-06,lwd=0.8)

text(1.5, 3.77e-06,"*",cex=2,col=("black"))
text(2, 4.77e-06, "*",cex=2,col=("black"))
text(5.5, 3.47e-06,"*",cex=2,col=("black"))
text(6, 4.47e-06, "*",cex=2,col=("black"))

legend("topleft", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
legend("bottomright", c("slow-med p* < 2.2e-16","slow-fast  p* < 2.2e-16","med-fast  p* < 2.2e-16 (TKO only)"), inset = .02, cex=.7)
dev.off()	

# --------------
# Violin plots | 
# --------------
#install.packages("vioplot")
library("vioplot")

# Pull
number <- nrow(merged_non05)
png(file = paste0(output, "5_Plots_without05/Violinplots_Pull2.png"))
vioplot(merged_non05$WT_Pull, merged_non05$TK0_Pull,
        h=0.3,
        drawRect=T,
        col = c("blue","red"),
        at = c(1,2),
        names = c("WT", "TKO"),
        ylab="Elongation rate (kb/min)",
        main = paste0("WT and H1-TKO elongation rates in ", number, " replicates"),
        border="black",
        plotCentre = "line",
        lineCol=NA,
        rectCol=NA,
        ylim=c(-2,12),
        boxwex = 0.2,
        areaEqual = F,
        xlab=paste("W = 13239520, p-value = 1.584e-11"),
        wex=0.5)
legend("bottomright", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
dev.off()

# replicates and quantiles
png(file = paste0(output, "5_Plots_without05/Violinplots_Pull_tertiles.png"))
vioplot(WTslow$WT_Pull, TKOslow$TK0_Pull, 
        WTmedium$WT_Pull, TKOmedium$TK0_Pull, 
        WTfast$WT_Pull, TKOfast$TK0_Pull,
        h=0.3,
        col = c("blue", "red"),
        at = c(1,2,4,5,7,8),
        names = c("","slow","med","","fast",""),
        ylab="Elongation rate (kb/min)",
        main = "WT and TKO elongation rates - replicates",
        border=1,
        rectCol=NA,
        lineCol=NA,
        plotCentre = "line",
        areaEqual = F)
legend("topleft", c("WT","TKO"), col=c("blue","red"), lwd=4, inset = .02)
legend("bottomright", c("slow p = 9.22e-10","med p* < 2.2e-16", "fast  p* < 2.2e-16"), inset = .02, cex=0.8)
dev.off()

# ----------------------
# Gene size histograms | 
# ----------------------
# Histogram All genes
data_by_size <- merged_non05[order(merged_non05$Size, decreasing = TRUE),] # sort by size
head(data_by_size)
data_size_mean <- round(mean(data_by_size[, 9]), digit=2) # get size mean
data_number <- nrow(data_by_size)
data_size_mean
data_number

png(file = paste0(output, "5_Plots_without05/Histograms_Pull_All_genes_size.png"))
hist(data_by_size[, 9],
     breaks=100,
     col=0,
     main = paste0("All genes (",data_number,") \n mean size = ", data_size_mean, " Kb"),
     xlab = paste("Gene size (kb)"))
dev.off()	

# Histogram WT slow genes
WTslow_by_size <- WTslow[order(WTslow$Size, decreasing = TRUE),] # sort by size
WTslow_size_mean <- round(mean(WTslow_by_size[, 9]), digit=2)   # get size mean
WTslow_number <- nrow(WTslow_by_size)

png(file = paste0(output, "5_Plots_without05/Histograms_Pull_WT_slow_genes_size_tertiles.png"))
hist(WTslow_by_size[, 9],
     breaks=100,
     col=4,
     main = paste0("Slow WT genes (",WTslow_number,") \n mean size = ", WTslow_size_mean, " Kb"),
     xlab = paste("Gene size (kb)"))
dev.off()	

# Histogram TKO slow genes
TKOslow_by_size <- TKOslow[order(TKOslow$Size, decreasing = TRUE),] # sort by size
TKOslow_size_mean <- round(mean(TKOslow_by_size[, 9]), digit=2)    # get size mean
TKOslow_size_mean
TKOslow_number <- nrow(TKOslow_by_size)

png(file = paste0(output, "5_Plots_without05/Histograms_Pull_TKO_slow_genes_size_tertiles.png"))
hist(TKOslow_by_size[, 9],
     breaks=100,
     col=2,
     main = paste0("Slow TKO genes (",TKOslow_number,") \n mean size = ", TKOslow_size_mean, " Kb"),
     xlab = paste("Gene size (kb)"))
dev.off()