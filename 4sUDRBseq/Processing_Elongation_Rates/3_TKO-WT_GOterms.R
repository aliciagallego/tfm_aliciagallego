#!/usr/bin/env Rscript

####################
## GO-terms graphs #
####################

# -------
# Paths |
# -------
GOterms_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/8_TKOvsWT_rates/Go-terms/"
output <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/8_TKOvsWT_rates/Go-terms_2SD/")

# -----------
# Open data |
# -----------
GOterms_list = list.files(GOterms_path, pattern="*2.txt")

for (i in seq_along(GOterms_list)) {
  filename <- sub("2.txt.*", "", GOterms_list[i])
  df <- read.table(paste0(GOterms_path,GOterms_list[i]), header = TRUE, sep="\t", stringsAsFactors=FALSE)
  #df <- df[-(1:11),] # remove first 11 rows
  df$FDR_log <- -log10(df[,8])
  df <- df[!(df[,5]=="-"),]
  df <- df[order(df$FDR_log, decreasing = T),]
  df[,1] <- sub(" \\(GO.*","", df[,1])
  assign(filename, df)

}

# -----------
# Bar plots |
# -----------

# All genes
png(file = paste0(output, "GO_AllGenes_biological_processes.png"))
ordered <- AllGenes_biological_processes[order(AllGenes_biological_processes[,9], decreasing = T),]
par(mar=c(5, 17, 3, 2))
barplot(ordered[,9][30:1],  
        space=c(0.5,0.5),
        legend.text=TRUE,
        beside=TRUE,
        horiz=TRUE,
        axes=TRUE,
        xlim=c(0,100),
        col="#69b3a2",
        border=NA,
        width=0.1,
        names.arg=ordered[,1][30:1], 
        cex.names=0.8, las=1,
        main="All genes - Biological processes",
        xlab = "-log10(FDR)")
dev.off()

png(file = paste0(output, "GO_AllGenes_molecular_function.png"))
ordered <- AllGenes_molecular_function[order(AllGenes_molecular_function[,9], decreasing = T),]
par(mar=c(5, 20, 3, 2))
barplot(ordered[,9][30:1],  
        space=c(0.5,0.5),
        legend.text=TRUE,
        beside=TRUE,
        horiz=TRUE,
        axes=TRUE,
        xlim=c(0,100),
        col="coral2",
        border=NA,
        width=0.1,
        names.arg=ordered[,1][30:1], 
        cex.names=0.8, las=1,
        main="All genes - Molecular functions",
        xlab = "-log10(FDR)")
dev.off()

# TKO Fast genes
png(file = paste0(output, "GO_TKOfast_biological_processes.png"))
ordered <- TKOfast_biological_processes[order(TKOfast_biological_processes[,9], decreasing = F),]
par(mar=c(5, 17, 3, 2))
barplot(ordered[,9][1:10], #FDR -log10 value 
        space=c(0.5,0.5),
        legend.text=TRUE,
        beside=TRUE,
        horiz=TRUE,
        axes=TRUE,
        xlim=c(0,4),
        col="#69b3a2",
        border=NA,
        width=0.1,
        names.arg=ordered[,1][1:10], 
        cex.names=0.8, las=1,
        main="TKO Fast - Biological processes",
        xlab = "-log10(FDR)")
#text("TKO Fast - Biological processes", x=3, y = 0.6, font=2, cex=1.2)
dev.off()

png(file = paste0(output, "GO_TKOfast_molecular_function.png"))
ordered <- TKOfast_molecular_function[order(TKOfast_molecular_function[,9], decreasing = F),]
head(TKOfast_molecular_function)
par(mar=c(5, 15, 3, 2))
barplot(ordered[,9],  
        space=c(0.5,0.5),
        legend.text=TRUE,
        beside=TRUE,
        horiz=TRUE,
        axes=TRUE,
        xlim=c(0,3.5),
        col="coral2",
        border=NA,
        width=0.1,
        names.arg=ordered[,1], 
        cex.names=0.8, las=1,
        main="TKO Fast - Molecular functions",
        xlab = "-log10(FDR)")
dev.off()

# TKO Fast+Slow genes
png(file = paste0(output, "GO_TKOfast_slow_biological_processes.png"))
ordered <- TKOfast_slow_biological_processes[order(TKOfast_slow_biological_processes[,9], decreasing = F),]
head(TKOfast_slow_biological_processes)
tail(TKOfast_slow_biological_processes)

par(mar=c(5, 17, 3, 2))
barplot(ordered[,9][1:19], #FDR -log10 value 
        space=c(0.5,0.5),
        legend.text=TRUE,
        beside=TRUE,
        horiz=TRUE,
        axes=TRUE,
        xlim=c(0,4),
        col="#69b3a2",
        border=NA,
        width=0.1,
        names.arg=ordered[,1][1:19], 
        cex.names=0.8, las=1,
        main="TKO Fast & Slow \nBiological processes",
        xlab = "-log10(FDR)")
#text("TKO Fast - Biological processes", x=3, y = 0.6, font=2, cex=1.2)
dev.off()

png(file = paste0(output, "GO_TKOfast_slow_molecular_function.png"))
ordered <- TKOfast_slow_molecular_function[order(TKOfast_slow_molecular_function[,9], decreasing = F),]
par(mar=c(5, 20, 3, 2))
barplot(ordered[,9][1:39],  
        space=c(0.5,0.5),
        legend.text=TRUE,
        beside=TRUE,
        horiz=TRUE,
        axes=TRUE,
        xlim=c(0,5.5),
        col="coral2",
        border=NA,
        width=0.1,
        names.arg=ordered[,1][1:39], 
        cex.names=0.8, las=1,
        main="TKO Fast & Slow \nMolecular functions",
        xlab = "-log10(FDR)")
dev.off()