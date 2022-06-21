#!/usr/bin/env Rscript

############################################################################
## Correlations elongation rates (4sUDRBseq or TTseq) between WT vs H1-TKO #
############################################################################

# -----------
# Libraries |
# -----------
library("ggpubr")
library("ggplot2")
library("reshape2")

# -------
# Paths |
# -------
TTseq_path <- "/path/4sUDRB/Elongation_rate/Rate_calculation/Elongation_rate_20Kb_Pull_processed_without05.txt"
TTseq_replicates_path <- "/path/4sUDRB/Elongation_rate/Rate_calculation/Elongation_rate_20Kb_Repl_processed_without05.txt"

output <- "/path/4sUDRB/Elongation_rate/Correlations/"
output2 <- "/path/4sUDRB/Elongation_rate/TKOvsWT_rates/"

# -----------
# Open data |
# -----------
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))

TTseq_replicates <- read.table(TTseq_replicates_path,h=T,sep="\t",stringsAsFactors=F)

# ------------------
# Correlation test |
# ------------------

# Pull
png(file = paste0(output, "TTseq_WT_TKO/Spearman_TTseq_pull_WT_TKO.png"))
TTseq$Gene_name <- rownames(TTseq)
ggscatter(TTseq, y = "TTseq_WT_pull", x = "TTseq_TKO_pull",
          color="darkmagenta",shape = 10,
          size=0.8,
          add.params = list(color = "darkmagenta", fill = "lightgray"),
          add = "reg.line", 
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "spearman",
          alpha=.3,
          ylim=c(0,12),
          xlim=c(0,12),
          xlab = "WT (Kb/min)", ylab = "TKO (Kb/min)") +
  theme_minimal()+
  geom_abline(intercept = 0, slope = 1, lwd=0.3, color="grey2")+
  stat_cor(method = "spearman", color = "darkmagenta", label.x = 7, label.y = 0.5,cex=5)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n elongation rates in WT vs TKO (", number, " genes)"))
dev.off()

# Replicates
png(file = paste0(output, "TTseq_WT_TKO/Spearman_TTseq_WT1_WT2.png"))
TTseq_replicatesWT <- TTseq_replicates[which(TTseq_replicates$WT1 > 0.5 & TTseq_replicates$WT2 > 0.5),]
number <- nrow(TTseq_replicatesWT)
ggscatter(TTseq_replicatesWT, x = "WT1", y = "WT2",
          color="blue",shape = 10,
          size=0.8,
          add.params = list(color = "blue", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = F, cor.method = "spearman",
          alpha=.3,
          ylim=c(0,12),
          xlim=c(0,12),
          xlab = "WT1 (Kb/min)", ylab = "WT2 (Kb/min)") +
  theme_minimal()+
  geom_abline(intercept = 0, slope = 1, lwd=0.3, color="grey2")+
  stat_cor(method = "spearman", color = "blue", label.x = 1, label.y = 11, cex=4)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n elongation rates in WT1 vs WT2 (", number, " genes)"))
dev.off()

png(file = paste0(output, "TTseq_WT_TKO/Spearman_TTseq_TKO1_TKO2.png"))
TTseq_replicatesTKO <- TTseq_replicates[which(TTseq_replicates$TKO1 > 0.5 & TTseq_replicates$TKO2 > 0.5),]
number <- nrow(TTseq_replicatesTKO)
ggscatter(TTseq_replicatesTKO, x = "TKO1", y = "TKO2",
          color="red",shape = 10,
          size=0.8,
          add.params = list(color = "red", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = F, cor.method = "spearman",
          alpha=.3,
          ylim=c(0,12),
          xlim=c(0,12),
          xlab = "TKO1 (Kb/min)", ylab = "TKO2 (Kb/min)") +
  theme_minimal()+
  geom_abline(intercept = 0, slope = 1, lwd=0.3, color="grey2")+
  stat_cor(method = "spearman", color = "red", label.x = 1, label.y = 11, cex=4)+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle(paste0("Spearman correlation \n elongation rates in TKO1 vs TKO2 (", number, " genes)"))
dev.off()
