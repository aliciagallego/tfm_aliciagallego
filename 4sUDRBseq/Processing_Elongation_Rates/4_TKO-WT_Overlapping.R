#!/usr/bin/env Rscript

###########################
## Overlap among tertiles # 
###########################

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull_processed_without05.txt"
output <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/8_TKOvsWT_rates/Overlaps/")

# -----------
# Open data |
# -----------
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT_pull","TTseq_TKO_pull","Size"))

# --------------------
# Add tertile labels |
# --------------------
QWT <- quantile(TTseq$TTseq_WT_pull, probs = seq(0, 1, 1/3), na.rm = FALSE)  
upWT <-  QWT[[3]]  
lowWT<- QWT[[2]]

QTKO <- quantile(TTseq$TTseq_TKO_pull, probs = seq(0, 1, 1/3), na.rm = FALSE)  
upTKO <-  QTKO[[3]]  
lowTKO<- QTKO[[2]]

TTseq$TTseq_WT_pull <lowWT

count <- 0

for(i in 1:nrow(TTseq)) {
  if(TTseq[i,7] <lowWT) {
    TTseq[i,10]="1WT_slow"
  }
  else if(TTseq[i,7] >= lowWT & TTseq[i,7] <=upWT) {
    TTseq[i,10]="2WT_medium"
  }
  else if (TTseq[i,7] >upWT){
    TTseq[i,10]="3WT_fast"
  }
  count = count+1
}
print(count)

count <- 0
for(i in 1:nrow(TTseq)) {
  if(TTseq[i,8] <lowTKO) {
    TTseq[i,11]="1TKO_slow"
  }
  else if(TTseq[i,8] >= lowTKO & TTseq[i,8] <=upTKO) {
    TTseq[i,11]="2TKO_medium"
  }
  else if (TTseq[i,8] >upTKO){
    TTseq[i,11]="3TKO_fast"
  }
  count = count+1
}
print(count)

colnames(TTseq)[10] <- "WT_label"
colnames(TTseq)[11] <- "TKO_label"

length(which(TTseq$WT_label == "1WT_slow"))
length(which(TTseq$WT_label == "2WT_medium"))
length(which(TTseq$WT_label == "3WT_fast"))

length(which(TTseq$TKO_label == "1TKO_slow"))
length(which(TTseq$TKO_label == "2TKO_medium"))
length(which(TTseq$TKO_label == "3TKO_fast"))

# ------------------------------------
# List of genes that change the most |
# ------------------------------------

WTfast_TKOslow <- TTseq[which(TTseq$WT_label == "3WT_fast" & TTseq$TKO_label == "1TKO_slow"),]
WTslow_TKOfast <- TTseq[which(TTseq$WT_label == "1WT_slow" & TTseq$TKO_label == "3TKO_fast"),]

WTfast_TKOslow <- WTfast_TKOslow[1]
WTslow_TKOfast <- WTslow_TKOfast[1]

# ------------------
# Save merged file |
# ------------------
write.table(WTfast_TKOslow, file = paste0(output,"TTseq_GeneList_WTfast_TKOslow.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

write.table(WTslow_TKOfast, file = paste0(output,"TTseq_GeneList_WTslow_TKOfast.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)

# -----------
# Bar plots |
# -----------
png(file = paste0(output, "Overlapping_TKO_across_WT.png"))
TKO <- ggplot(data = TTseq) + 
  geom_bar(mapping = aes(x = WT_label, fill = TKO_label),width = 0.6)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="right")+
  ggtitle(paste0("TKO genes across WT groups"))+
  xlab(element_blank())+ylab(element_blank())
TKO + scale_fill_manual(values=c("firebrick1", "firebrick3", "firebrick4"))
dev.off()

png(file = paste0(output, "Overlapping_WT_across_TKO.png"))
WT <- ggplot(data = TTseq) + 
  geom_bar(mapping = aes(x = TKO_label, fill = WT_label),width = 0.6)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position="right")+
  ggtitle(paste0("WT genes across TKO groups"))+
  xlab(element_blank())+ylab(element_blank())
WT + scale_fill_manual(values=c("steelblue1", "steelblue3", "steelblue4"))
dev.off()