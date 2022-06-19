#!/usr/bin/env Rscript

############################################################
## Correlations elongation rates Maslon and m6A (meRIPseq) # 20220122
############################################################

# -----------
# Libraries |
# -----------
#library("ggpubr")

# -------
# Paths |
# -------
maslon_m6A_path <- "/media/cc/A/Alicia/Maslon2019/Maslon_output/maslon_rates_m6A.txt"
output <- "/media/cc/A/Alicia/NGS/Correlations/elongationMaslon_m6A/"

# -----------
# Open data |
# -----------
maslon_m6A <- read.table(maslon_m6A_path,h=T,sep="\t",stringsAsFactors=FALSE)
nrow(maslon_m6A)

# ---------------------------------------------------------------
# Remove outliers from meRIP data and log transform RNApol data |
# ---------------------------------------------------------------
Q <- quantile(maslon_m6A$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(maslon_m6A$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(maslon_m6A, maslon_m6A$meRIP_WT > low & maslon_m6A$meRIP_WT < up)
head(saved_values_WT)
nrow(saved_values_WT)

# ------------------
# Correlation test |
# ------------------
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson", "kendall", "spearman"))
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson"))

# Pearson resulted to be more significant
png(file = paste0(output, "Pearson_MaslonElongation_meRIP_WT.png"))
ggscatter(saved_values_WT, y = "meRIP_WT", x = "elongation_b_min",
          color="lightblue3",shape = 20,
          size=1.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "m6A levels", xlab = "Elongation rate Maslon data (Kb/min)") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Pearson correlation m6A levels vs. elongation rates (Maslon data) - WT")
dev.off()