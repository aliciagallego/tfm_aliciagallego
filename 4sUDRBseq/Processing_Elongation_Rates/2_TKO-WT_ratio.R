#!/usr/bin/env Rscript

##########################################################
## Ratio elongation rates (4sUDRBseq or TTseq) TKO vs WT #
##########################################################

# -----------
# Libraries |
# -----------
library("ggpubr")
library("ggplot2")
library("reshape2")
library("RColorBrewer")

# -------
# Paths |
# -------
TTseq_path <- "/path/4sUDRB/Elongation_rate/Rate_calculation/Elongation_rate_20Kb_Pull_processed_without05.txt"

output2 <- "/path/4sUDRB/Elongation_rate/Correlations/"
output3 <- "/path/4sUDRB/Elongation_rate/TKOvsWT_rates/TKO_comparisons/"

# -----------
# Open data |
# -----------
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))

# -----------------
# Ratios TKO / WT |
# -----------------
TTseq$TKO_WT_ratio <- TTseq$TTseq_TKO_pull/TTseq$TTseq_WT_pull
TTseq$log_ratio <- log(TTseq$TKO_WT_ratio)
TTseq <- TTseq[order(TTseq$log_ratio, decreasing = TRUE),]

# -----------------------------
# Mean and Standard deviation |
# -----------------------------
mean <- mean(TTseq$log_ratio)
sd <- sd(TTseq$log_ratio)

# ----------
# Q-Q plot | confirm ratios have a normal distribution
# ----------
ggqqplot(TTseq$log_ratio)

# ---------------------------------------
# TKO genes significantly faster/slower |
# ---------------------------------------
TTseq_gene_list <- TTseq[1]

TTseq_candidates_all <- TTseq[which(TTseq$log_ratio >= (mean+2*sd) | TTseq$log_ratio <= (mean-2*sd)),] #2SD (95% confianza)
TTseq_candidates_all_1SD <- TTseq[which(TTseq$log_ratio >= (mean+sd) | TTseq$log_ratio <= (mean-sd)),] #1SD (68% confianza)
TTseq_candidates_all_1SD <- TTseq_candidates_all_1SD[1]

TTseq_candidates_fast <- TTseq[which(TTseq$log_ratio >= (mean+2*sd)),] #+2SD (95% confianza)
TTseq_candidates_fast <- TTseq_candidates_fast[1]

TTseq_candidates_fast_1SD <- TTseq[which(TTseq$log_ratio >= (mean+sd)),] #+1SD (68% confianza)
TTseq_candidates_fast_1SD <- TTseq_candidates_fast_1SD[1]

TTseq_candidates_slow <- TTseq[which(TTseq$log_ratio <= (mean-2*sd)),] #-2SD (95% confianza)
TTseq_candidates_slow <- TTseq_candidates_slow[1]

TTseq_candidates_slow_1SD <- TTseq[which(TTseq$log_ratio <= (mean-sd)),] #-1SD (68% confianza)
TTseq_candidates_slow_1SD <- TTseq_candidates_slow_1SD[1]

# ------------------
# Save merged file |
# ------------------
write.table(TTseq_gene_list, file = paste0(output2,"TTseq_GeneList.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(TTseq_candidates_all, file = paste0(output2,"TTseq_TKO_fast&slow.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(TTseq_candidates_fast, file = paste0(output2,"TTseq_TKO_fast_GeneList.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(TTseq_candidates_slow, file = paste0(output2,"TTseq_TKO_slow_GeneList.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(TTseq_candidates_all_1SD, file = paste0(output2,"TTseq_TKO_fast_slow_1SD.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(TTseq_candidates_fast_1SD, file = paste0(output2,"TTseq_TKO_fast_GeneList_1SD.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(TTseq_candidates_slow_1SD, file = paste0(output2,"TTseq_TKO_slow_GeneList_1SD.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# ---------------------------------------------------
# TKO genes significantly faster/slower with labels |
# ---------------------------------------------------
count <- 0
for(i in 1:nrow(TTseq)) {
  if(TTseq[i,11] >= (mean+2*sd) | TTseq[i,11] <= (mean-2*sd)) {
    TTseq[i,12]="±2SD"
  }
  else if((TTseq[i,11] >= (mean+sd) & TTseq[i,11] < (mean+2*sd)) | (TTseq[i,11] <= (mean-sd) & TTseq[i,11] > (mean-2*sd))) {
    TTseq[i,12]="±1SD"
  }
  else {
    TTseq[i,12]="68%"
  }
  count = count+1
}
print(count)

colnames(TTseq)[12] <- "CI"
length(which(TTseq$CI == "±2SD"))
findTermsInFileNames(TTseq, "±2SD")

# ------------
# Histograms |
# ------------
hist(TTseq$TKO_WT_ratio,
     breaks=50,
     col="aquamarine3",
     border="aquamarine3",
     main = paste0("Ratio TKO/WT rates"),
     xlab = paste("ratio \n bin_size = 50"))
abline(v = 1, col="black", lwd=1, lty=2)

hist(TTseq$log_ratio,
     breaks=150,
     col="aquamarine4",
     border="aquamarine4",
     main = paste0("log Ratio TKO/WT rates"),
     xlab = paste("log ratio \n bin_size = 150"))
abline(v = 1, col="black", lwd=1, lty=2)

# ---------------
# Density plots | 
# ---------------
png(file = paste0(output2, "Density_plot_TKO_WT_ratio_log.png"))
plot(density(TTseq$log_ratio, bw=.05), 
     col="aquamarine4",  
     main ="Log Elongation rate TKO/WT ratio", 
     cex.main=1, 
     lwd=2,
     xlab = "log TKO/WT\nbw = 0.05")
abline(v=c(-2*sd,2*sd),col="black", lwd=1, lty=2)
legend("topright", paste0("CI: ± 2SD"), 
       inset = .02,cex=0.9)
dev.off()

number<-nrow(TTseq_candidates)
plot(density(TTseq_candidates$log_ratio,bw=.2),
     col="salmon",
     lwd=2,
     main = paste0("TKO genes with faster/slower significant rates\n(", number, " genes, threshold = 95%)"),
     xlab="log TKO/WT ratio")

# --------------
# Scatter plot | 
# --------------
png(file = paste0(output2, "Scatter_plot_TKO_WT_ConfidenceInterval.png"))
number <- nrow(TTseq)
ggscatter(TTseq, y = "TTseq_WT_pull", x = "TTseq_TKO_pull",
          color="CI",
          alpha=.5,
          shape = 10,
          palette = c("#008080","#DC143C","black"),
          size=0.8,
          ylim=c(0.5,12),
          xlim=c(0.5,12),
          xlab = "WT (Kb/min)", ylab = "TKO (Kb/min)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(paste0("Elongation rates in WT and TKO\n(",number," genes)"))+
  geom_abline(intercept = 0, slope = c(1.65,2.9,0.35,0.58), color = c("#008B8B","red","red","#008B8B"),linetype=2,lwd=1)
dev.off()

# ------------------------------
# Fast/Slow TKO Genes and size |
# ------------------------------
par(mar=c(5,5,5,5))

png(file = paste0(output3, "Size_fast.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_ratio>= (mean+2*sd)),]$Size, TTseq[which(TTseq$log_ratio <(mean+2*sd)),]$Size)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_ratio>= (mean+2*sd)),]$Size),3)
mean2<-round(mean(TTseq[which(TTseq$log_ratio <(mean+2*sd)),]$Size),3)
boxplot(TTseq[which(TTseq$log_ratio>= (mean+2*sd)),]$Size, TTseq[which(TTseq$log_ratio <(mean+2*sd)),]$Size,
        col = c("darkgoldenrod1","darkgoldenrod3"),
        at = c(1,2),
        names = c("TKO fast","Rest"),
        ylab="Gene size (Kb)",
        main = paste0("Gene size of TKO faster genes"),
        xlab=paste("W = 556988, p-val = 0.003"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "Size_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_ratio<= (mean-2*sd)),]$Size, TTseq[which(TTseq$log_ratio >(mean-2*sd)),]$Size)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_ratio<= (mean-2*sd)),]$Size),3)
mean2<-round(mean(TTseq[which(TTseq$log_ratio >(mean-2*sd)),]$Size),3)
boxplot(TTseq[which(TTseq$log_ratio<= (mean-2*sd)),]$Size, TTseq[which(TTseq$log_ratio >(mean-2*sd)),]$Size,
        col = c("lightcyan3","lightcyan4"),
        at = c(1,2),
        names = c("TKO slow","Rest"),
        ylab="Gene size (Kb)",
        main = paste0("Gene size of TKO slower genes"),
        xlab=paste("W = 722184, p-val = 6.8e-05"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "Size_fast_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_ratio>= (mean+2*sd) | TTseq$log_ratio<= (mean-2*sd)),]$Size, TTseq[which(TTseq$log_ratio <(mean+2*sd)|TTseq$log_ratio >(mean-2*sd)),]$Size)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_ratio>= (mean+2*sd) | TTseq$log_ratio<= (mean-2*sd)),]$Size),3)
mean2<-round(mean(TTseq[which(TTseq$log_ratio <(mean+2*sd)|TTseq$log_ratio >(mean-2*sd)),]$Size),3)
boxplot(TTseq[which(TTseq$log_ratio>= (mean+2*sd) | TTseq$log_ratio<= (mean-2*sd)),]$Size, TTseq[which(TTseq$log_ratio <(mean+2*sd)|TTseq$log_ratio >(mean-2*sd)),]$Size,
        col = c("springgreen3","springgreen4"),
        at = c(1,2),
        names = c("TKO fast & slow","Rest"),
        ylab="Gene size (Kb)",
        main = paste0("Gene size of TKO faster and slower genes"),
        xlab=paste("W = 1327862, p-value = 3.5e-06"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

# -------------------------------------
# Fast/Slow TKO Genes and exon number |
# -------------------------------------
par(mar=c(5,5,5,5))
head(TTseq)
png(file = paste0(output3, "Exon_fast.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_ratio>= (mean+2*sd)),]$Exon_number, TTseq[which(TTseq$log_ratio <(mean+2*sd)),]$Exon_number)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_ratio>= (mean+2*sd)),]$Exon_number),3)
mean2<-round(mean(TTseq[which(TTseq$log_ratio <(mean+2*sd)),]$Exon_number),3)
boxplot(TTseq[which(TTseq$log_ratio>= (mean+2*sd)),]$Exon_number, TTseq[which(TTseq$log_ratio <(mean+2*sd)),]$Exon_number,
        col = c("darkgoldenrod1","darkgoldenrod3"),
        at = c(1,2),
        names = c("TKO fast","Rest"),
        ylab="Exon number",
        main = paste0("Exon number of TKO faster genes"),
        xlab=paste("W = 483472, p-val = 0.5745"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "Exon_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_ratio<= (mean-2*sd)),]$Exon_number, TTseq[which(TTseq$log_ratio >(mean-2*sd)),]$Exon_number)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_ratio<= (mean-2*sd)),]$Exon_number),3)
mean2<-round(mean(TTseq[which(TTseq$log_ratio >(mean-2*sd)),]$Exon_number),3)
boxplot(TTseq[which(TTseq$log_ratio<= (mean-2*sd)),]$Exon_number, TTseq[which(TTseq$log_ratio >(mean-2*sd)),]$Exon_number,
        col = c("lightcyan3","lightcyan4"),
        at = c(1,2),
        names = c("TKO slow","Rest"),
        ylab="W = 616478, p-val = 0.6287",
        main = paste0("Exon number of TKO slower genes"),
        xlab=paste("W = 616478, p-val = 0.6287"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	

png(file = paste0(output3, "Exon_fast_slow.png"))
wilcox <- wilcox.test(TTseq[which(TTseq$log_ratio>= (mean+2*sd) | TTseq$log_ratio<= (mean-2*sd)),]$Exon_number, TTseq[which(TTseq$log_ratio <(mean+2*sd)|TTseq$log_ratio >(mean-2*sd)),]$Exon_number)
wilcox
mean1<-round(mean(TTseq[which(TTseq$log_ratio>= (mean+2*sd) | TTseq$log_ratio<= (mean-2*sd)),]$Exon_number),3)
mean2<-round(mean(TTseq[which(TTseq$log_ratio <(mean+2*sd)|TTseq$log_ratio >(mean-2*sd)),]$Exon_number),3)
boxplot(TTseq[which(TTseq$log_ratio>= (mean+2*sd) | TTseq$log_ratio<= (mean-2*sd)),]$Exon_number, TTseq[which(TTseq$log_ratio <(mean+2*sd)|TTseq$log_ratio >(mean-2*sd)),]$Exon_number,
        col = c("springgreen3","springgreen4"),
        at = c(1,2),
        names = c("TKO fast & slow","Rest"),
        ylab="Exon number",
        main = paste0("Exon number of TKO faster and slower genes"),
        xlab=paste("W = 1148641, p-val = 0.4894"),
        outline = F,
        boxwex = 0.3)
text(mean1, x=1, y=mean1, col= "black")
text(mean2, x=2, y=mean2, col= "black")
dev.off()	
