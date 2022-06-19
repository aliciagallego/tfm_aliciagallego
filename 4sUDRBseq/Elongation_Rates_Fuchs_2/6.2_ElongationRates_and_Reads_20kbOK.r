#!/usr/bin/env Rscript

########################################################
## Get number of reads and elongation rates from TTseq # 20220328
########################################################

# -------
# Paths |
# -------
TTseq_Intersect_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/2_Intersect_RefSeqBED20Kb_BAM/"
TTseq_Rates_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/1.1_Rate_calculation/Elongation_rate_5min_20220401_20Kb_size.txt"
output <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/"

# -----------
# Open data |
# -----------
TTseq_Rates <- read.table(TTseq_Rates_path,h=T,sep="\t",stringsAsFactors=FALSE)
nrow(TTseq_Rates)

TTseq_list = list.files(TTseq_Intersect_path, pattern="*.bed")
TTseq_list

for (i in seq_along(TTseq_list)) {
  filename <- sub(".bed", "", TTseq_list[i])
  filename <- sub("-", "_", filename)
  df <- read.table(paste0(TTseq_Intersect_path,TTseq_list[i]), header = FALSE, sep="\t", stringsAsFactors=FALSE,
                   col.names = c("Chr","Start","End","Gene_name", "NA1", "Strand",paste0("TTseq_",filename)))
  df <- subset(df, select=-c(NA1)) # remove uninformative columns
  assign(filename, df)
}
nrow(MG9_12)
head(MG9_12)

# --------------------------------------------------
# Remove white spaces in Gene_names in TTseq Rates | 
# --------------------------------------------------
TTseq_Rates$Gene_name <- gsub("\\s", "", TTseq_Rates$Gene_name)
TTseq_Rates[ TTseq_Rates == "NaN" ] <- NA

# ------------------------------
# Rates from bp/5min to Kb/min |
# ------------------------------
TTseq_Rates[,2] = TTseq_Rates[,2]/5000
TTseq_Rates[,3] = TTseq_Rates[,3]/5000
TTseq_Rates[,4] = TTseq_Rates[,4]/5000  
TTseq_Rates[,5] = TTseq_Rates[,5]/5000  
head(TTseq_Rates)

# -------------------
# Delete NaN values |
# -------------------
TTseq_Rates_noNA <- na.omit(TTseq_Rates)  
nrow(TTseq_Rates_noNA)

# ------------
# Read rates |
# ------------
TTseq_Rates_noNA$WTmean = (TTseq_Rates_noNA$WT1+TTseq_Rates_noNA$WT2)/2
TTseq_Rates_noNA$TKOmean = (TTseq_Rates_noNA$TKO1+TTseq_Rates_noNA$TKO2)/2
head(TTseq_Rates_noNA)

# ---------------------------------------------------
# Merge Elongation rates list and Input RefSeq list |
# ---------------------------------------------------
merged_Intersect<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                list(MG9_12,
                     MG9_14[,c("TTseq_MG9_14","Gene_name")],
                     MG9_16[,c("TTseq_MG9_16","Gene_name")],
                     MG9_18[,c("TTseq_MG9_18","Gene_name")]))

head(merged_Intersect)
nrow(merged_Intersect)

merged<- Reduce(function(x,y) merge(x = x, y = y, by="Gene_name", sort = F),
                list(merged_Intersect,
                     TTseq_Rates_noNA[,c("Gene_name","WT1", "WT2", "TKO1", "TKO2","WTmean","TKOmean")]))

head(merged)
nrow(merged)
# ------------
# Mean rates |
# ------------
merged$WT_reads_mean = (merged$TTseq_MG9_12+merged$TTseq_MG9_14)/2
merged$TKO_reads_mean = (merged$TTseq_MG9_16+merged$TTseq_MG9_18)/2
head(merged)

# -------------------------------------------
# Plot elongation rates and number of reads |
# -------------------------------------------
number <- nrow(merged)

png(file = paste0(output, "5_Plots/Plots_20Kb_list_withSize_ratios/Plot_rates_and_reads_WT.png"))
merged_ordered_WTrate <- merged[order(merged$WTmean,decreasing = F),] 
head(merged_ordered_WTrate)
plot(merged_ordered_WTrate$WT_reads_mean/1000,
     main=paste0("Elongation rates and reads in WT \n (", number, " genes)"),
     xlab="# genes",
     ylab="elongation rates & number of reads",
     pch=1, 
     cex=.6, 
     ylim = c(0,11),
     col="lightblue")
lines(merged_ordered_WTrate$WTmean,
      col="blue",
      lwd = 2,
      type = "l")
legend("topleft", inset=.02, legend=c("Number of reads/1000","Elongation rates Kb/min"),
       col=c("lightblue","blue"), lty=c(NA, 1),pch=c(1, NA), lwd=c(1,2), cex=0.8)
dev.off()	

png(file = paste0(output, "5_Plots/Plots_20Kb_list_withSize_ratios/Plot_rates_and_reads_TKO.png"))
merged_ordered_TKOrate <- merged[order(merged$TKOmean,decreasing = F),] 
head(merged_ordered_TKOrate)
plot(merged_ordered_TKOrate$TKO_reads_mean/1000,
     main=paste0("Elongation rates and reads in TKO \n (", number, " genes)"),
     xlab="# genes",
     ylab="elongation rates & number of reads",
     pch=1, 
     cex=.6, 
     ylim = c(0,11),
     col="salmon")
lines(merged_ordered_TKOrate$TKOmean,
     col="red",
     lwd = 2,
     type = "l")
legend("topleft", inset=.02, legend=c("Number of reads/1000","Elongation rates Kb/min"),
       col=c("salmon","red"), lty=c(NA, 1),pch=c(1, NA), lwd=c(1,2), cex=0.8)
dev.off()	

png(file = paste0(output, "5_Plots/Plots_20Kb_list_withSize_ratios/Plot_reads_and_rates_WT.png"))
merged_ordered_WTreads <- merged[order(merged$WT_reads_mean,decreasing = F),] 
head(merged_ordered_WTreads)
plot(merged_ordered_WTreads$WTmean,
     pch=16, 
     cex=.6, 
     col="lightblue",
     xlab="# genes", ylab="elongation rates & number of reads",
     #ylim=c(0,80),
     main="Elongation rates and reads in WT")
lines(merged_ordered_WTreads$WT_reads_mean/1000,
      col="blue",
      lwd = 2,
      type = "l")
legend("topleft", inset=.02, legend=c("Elongation rates Kb/min", "Number of reads/1000"),
       col=c("lightblue","blue"), lty=c(NA, 1),pch=c(16, NA), lwd=c(1,2), cex=0.8)
dev.off()	

png(file = paste0(output, "5_Plots/Plots_20Kb_list_withSize_ratios/Plot_reads_and_rates_TKO.png"))
merged_ordered_TKOreads <- merged[order(merged$TKO_reads_mean,decreasing = F),] 
head(merged_ordered_TKOreads)
plot(merged_ordered_TKOreads$WTmean,
     pch=16, 
     cex=.6, 
     col="lightsalmon",
     xlab="# genes", ylab="elongation rates & number of reads",
     main="Elongation rates and reads in TKO")
lines(merged_ordered_TKOreads$TKO_reads_mean/1000,
      col="red",
      lwd = 2,
      type = "l")
legend("topleft", inset=.02, legend=c("Elongation rates Kb/min", "Number of reads/1000"),
       col=c("lightsalmon","red"), lty=c(NA, 1),pch=c(16, NA), lwd=c(1,2), cex=0.8)
dev.off()	