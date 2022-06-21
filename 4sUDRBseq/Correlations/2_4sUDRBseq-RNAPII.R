#!/usr/bin/env Rscript

##########################################################################
## Correlations elongation rates (4sUDRBseq or TTseq) and RNAPII ChiPseq #
##########################################################################

# -----------
# Libraries |
# -----------
library("ggplot2")
library("ggpubr")
library("reshape2")

# -------
# Paths |
# -------
# RNAPII levels in GB and PR calculated in the LongList transcriptome (20621 genes) from RefSeq using Josemi's definition:
# GB: start coord + 500 bp / end coord with no changes
# PR: start coord - 500 bp / start coord + 500 bp

RNApolGB_path <- "/path/RNAPII/Normalized_data/GeneBody/RNAPII_RefSeq_LongList_GB_500bp.bed"
RNApolPR_path <- "/path/RNAPII/Normalized_data/Promoters/RNAPII_RefSeq_LongList_PR_500bp.bed"
TTseq_path <- "/path/4sUDRB/Elongation_rate/Rate_calculation/Elongation_rate_20Kb_Pull_processed_without05.txt"

output <- ("/path/4sUDRB/Elongation_rate/Correlations/")
output3 <- ("/path/4sUDRB/Elongation_rate/TKOvsWT_rates/TKO_comparisons/")

# -----------
# Open data |
# -----------
RNApolGB <- read.table(RNApolGB_path, h=T, sep="\t", stringsAsFactors=F)
RNApolPR <- read.table(RNApolPR_path, h=T, sep="\t", stringsAsFactors=F)
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))

# -------------
# Merge files |
# -------------
merged_no05<- Reduce(function(x,y) merge(x = x, y = y, by = "Gene_name", sort = F),
                list(RNApolGB[,c("Chr","Start","End","Strand","RNAPII_WTmean_GB","RNAPII_TKOmean_GB", "Gene_name")],
                     RNApolPR[,c("RNAPII_WTmean_PR","RNAPII_TKOmean_PR", "Gene_name")],
                     TTseq[,c("TTseq_WT_pull", "TTseq_TKO_pull", "Exon_number","Size","Gene_name")]))

# ------------------
# Save merged file |
# ------------------
write.table(merged_no05, file = paste0(output,"TTseq_RNApol/TTseq_RNApol_500bp.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# ------------------
# Correlation test |
# ------------------

# RNAPII GB WTmean - TTseq WT Pull > 0.5 values 
png(file = paste0(output, "TTseq_RNApol/Spearman_RNApolII_TTseq_WT_GB_no05.png"))
number_no05 <- nrow(merged_no05)
ggscatter(merged_no05, x = "TTseq_WT_pull", y = "RNAPII_WTmean_GB",
          color="blue",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "blue", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylim=c(0,3e-06),
          xlim=c(0,9),
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "blue", label.x = 5, label.y = 2.8e-06,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n RNApolII levels vs. elongation rates in WT GB \n (", number_no05, " genes)"))
dev.off()

# RNAPII GB TKOmean - TTseq TKO Pull > 0.5 values 
png(file = paste0(output, "TTseq_RNApol/Spearman_RNApolII_TTseq_TKO_GB_no05.png"))
number_no05 <- nrow(merged_no05)
ggscatter(merged_no05, x = "TTseq_TKO_pull", y = "RNAPII_TKOmean_GB",
          color="red2",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "red2", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylim=c(0,3e-06),
          xlim=c(0,9),
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "red2", label.x = 5, label.y = 2.8e-06,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n RNApolII levels vs. elongation rates in TKO GB \n (", number_no05, " genes)"))
dev.off()

# RNAPII PR WTmean - TTseq WT Pull > 0.5 values 
png(file = paste0(output, "TTseq_RNApol/Spearman_RNApolII_TTseq_WT_PR_no05.png"))
number <- nrow(merged_no05)
ggscatter(merged_no05, x = "TTseq_WT_pull", y = "RNAPII_WTmean_PR",
          color="blue",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "blue", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylim=c(0,1.5e-05),
          #xlim=c(0,9),
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "blue", label.x = 5.5, label.y = 1.25e-05,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n RNApolII levels vs. elongation rates in WT PR \n (", number, " genes)"))
dev.off()

# RNAPII PR TKOmean - TTseq TKO Pull > 0.5 values 
png(file = paste0(output, "TTseq_RNApol/Spearman_RNApolII_TTseq_TKO_PR_no05.png"))
number <- nrow(merged_no05)
ggscatter(merged_no05, x = "TTseq_TKO_pull", y = "RNAPII_TKOmean_PR",
          color="red2",shape = 10,
          size=0.5,
          alpha=.2,
          add.params = list(color = "red2", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylim=c(0,1.5e-05),
          #xlim=c(0,9),
          xlab = "Elongation rate (Kb/min)", ylab = "RNApolII") +
  theme_minimal()+
  stat_cor(method = "spearman", color = "red2", label.x = 5.5, label.y = 1.25e-05,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n RNApolII levels vs. elongation rates in TKO PR \n (", number, " genes)"))
dev.off()

# -------------------------------
# Ratio TKO/WT rates and RNAPII |
# -------------------------------
TTseq <- merged_no05
TTseq$rates_ratio <- TTseq$TTseq_TKO_pull/TTseq$TTseq_WT_pull
TTseq$RNAPIIGB_ratio <- TTseq$RNAPII_TKOmean_GB/TTseq$RNAPII_WTmean_GB
TTseq$RNAPIIPR_ratio <- TTseq$RNAPII_TKOmean_PR/TTseq$RNAPII_WTmean_PR
TTseq$log_rates <- log(TTseq$rates_ratio)
TTseq <- TTseq[order(TTseq$log_ratio, decreasing = TRUE),]

mean <- mean(TTseq$log_rates)
sd <- sd(TTseq$log_rates)

TTseq[is.na(TTseq)] <- 0 # In case 0/0 produced Nan values

for(i in 1:nrow(TTseq)) {
  if(TTseq[i,17] >= (mean+2*sd) | TTseq[i,17] <= (mean-2*sd)){
    TTseq[i,18]="±2SD"
  }
  else {
    TTseq[i,18]="95%"
  }
}

colnames(TTseq)[18] <- "CI"

# Gene body (GB)
png(file = paste0(output3, "RNAPIIGBratio_fast.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$RNAPIIGB_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)),]$RNAPIIGB_ratio)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$RNAPIIGB_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates <(mean+2*sd)),]$RNAPIIGB_ratio),3)
boxplot(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$RNAPIIGB_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)),]$RNAPIIGB_ratio,
        col = c("darkgoldenrod1","darkgoldenrod3"),
        at = c(1,2),
        names = c("TKO fast","Rest"),
        ylab="RNAPII GB TKO/WT ratio",
        main = paste0("RNAPII GB TKO/WT of TKO faster genes"),
        xlab=paste("W = 542725, p-val = 0.02333"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "RNAPIIGBratio_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$RNAPIIGB_ratio, TTseq[which(TTseq$log_rates >(mean-2*sd)),]$RNAPIIGB_ratio)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$RNAPIIGB_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates >(mean-2*sd)),]$RNAPIIGB_ratio),3)
boxplot(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$RNAPIIGB_ratio, TTseq[which(TTseq$log_rates >(mean-2*sd)),]$RNAPIIGB_ratio,
        col = c("lightcyan3","lightcyan4"),
        at = c(1,2),
        names = c("TKO slow","Rest"),
        ylab="RNAPII GB TKO/WT ratio",
        main = paste0("RNAPII GB TKO/WT ratio of TKO slower genes"),
        xlab=paste("W = 699535, p-val = 0.002383"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2+.01, col= "black")
dev.off()	

png(file = paste0(output3, "RNAPIIGBratio_fast_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$RNAPIIGB_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$RNAPIIGB_ratio)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$RNAPIIGB_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$RNAPIIGB_ratio),3)
boxplot(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$RNAPIIGB_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$RNAPIIGB_ratio,
        col = c("springgreen3","springgreen4"),
        at = c(1,2),
        names = c("TKO fast & slow","Rest"),
        ylab="RNAPII GB TKO/WT ratio",
        main = paste0("RNAPII GB TKO/WT ratio of TKO faster and slower genes"),
        xlab=paste("W = 1290950, p-val = 0.0003774"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2+.01, col= "black")
dev.off()	

# Promoter (PR)
head(TTseq)
png(file = paste0(output3, "RNAPIIPRratio_fast.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$RNAPIIPR_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)),]$RNAPIIPR_ratio)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$RNAPIIPR_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates <(mean+2*sd)),]$RNAPIIPR_ratio),3)
boxplot(TTseq[which(TTseq$log_rates>= (mean+2*sd)),]$RNAPIIPR_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)),]$RNAPIIPR_ratio,
        col = c("darkgoldenrod1","darkgoldenrod3"),
        at = c(1,2),
        names = c("TKO fast","Rest"),
        ylab="RNAPII PR TKO/WT ratio",
        main = paste0("RNAPII PR TKO/WT of TKO faster genes"),
        xlab=paste("W = 492894, p-val = 0.9175"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "RNAPIIPRratio_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$RNAPIIPR_ratio, TTseq[which(TTseq$log_rates >(mean-2*sd)),]$RNAPIIPR_ratio)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$RNAPIIPR_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates >(mean-2*sd)),]$RNAPIIPR_ratio),3)
boxplot(TTseq[which(TTseq$log_rates<= (mean-2*sd)),]$RNAPIIPR_ratio, TTseq[which(TTseq$log_rates >(mean-2*sd)),]$RNAPIIPR_ratio,
        col = c("lightcyan3","lightcyan4"),
        at = c(1,2),
        names = c("TKO slow","Rest"),
        ylab="RNAPII PR TKO/WT ratio",
        main = paste0("RNAPII PR TKO/WT ratio of TKO slower genes"),
        xlab=paste("W = 657021, p-value = 0.2147"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2+.01, col= "black")
dev.off()	

png(file = paste0(output3, "RNAPIIPRratio_fast_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$RNAPIIPR_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$RNAPIIPR_ratio)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$RNAPIIPR_ratio),3)
mean2<-round(mean(TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$RNAPIIPR_ratio),3)
boxplot(TTseq[which(TTseq$log_rates>= (mean+2*sd) | TTseq$log_rates<= (mean-2*sd)),]$RNAPIIPR_ratio, TTseq[which(TTseq$log_rates <(mean+2*sd)|TTseq$log_rates >(mean-2*sd)),]$RNAPIIPR_ratio,
        col = c("springgreen3","springgreen4"),
        at = c(1,2),
        names = c("TKO fast & slow","Rest"),
        ylab="RNAPII PR TKO/WT ratio",
        main = paste0("RNAPII PR TKO/WT ratio of TKO faster and slower genes"),
        xlab=paste("W = 1198605, p-value = 0.4189"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2+.01, col= "black")
dev.off()	
