#!/usr/bin/env Rscript

####################
## Metaplots Pulls #
####################

# -------
# Paths |
# -------
WT_pull_path <- "/path/4suDRB/Pull/Histograms/MG9-WT_Pull_matrix_sum_001.txt"
TKO_pull_path <-"/path/4suDRB/Pull/Histograms/MG9-TKO_Pull_matrix_sum_001.txt"

output <- ("path/4suDRB/Pull/Histograms/")

# -----------
# Open data |
# -----------
WT_pull <- read.table(WT_pull_path,h=T,sep="\t",stringsAsFactors=FALSE)
TKO_pull <- read.table(TKO_pull_path,h=T,sep="\t",stringsAsFactors=FALSE)

# ----------
# Metaplot |
# ----------

# WT Pull
png(file = paste0(output, "Metaplot_TTseq_WT_Pull.png"))
plot(WT_pull$X, WT_pull$WT5_mean,type="l", ylim=c(0,60000),xlim=c(-2000,50000),lwd=2, col="blue",
     xlab = "Gene body (-2Kb to +50Kb)",
     ylab = "Number of reads")
lines(WT_pull$X,WT_pull$WT0_mean, col="grey",lwd=2)
legend("topright", c("WT 0'","WT 5'"), 
       col=c("grey","blue"), lwd=2, inset =.02, cex=0.9)
dev.off()

# TKO Pull
png(file = paste0(output, "Metaplot_TTseq_TKO_Pull.png"))
plot(TKO_pull$X, TKO_pull$TKO5_mean,type="l", ylim=c(0,60000),xlim=c(-2000,50000),lwd=2, col="red",
     xlab = "Gene body (-2Kb to +50Kb)",
     ylab = "Number of reads")
lines(TKO_pull$X,TKO_pull$TKO0_mean, col="grey",lwd=2)
legend("topright", c("TKO 0'","TKO 5'"), 
       col=c("grey","red"), lwd=2, inset =.02, cex=0.9)
dev.off()
