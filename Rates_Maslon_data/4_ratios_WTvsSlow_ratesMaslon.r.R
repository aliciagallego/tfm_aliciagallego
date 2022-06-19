#!/usr/bin/env Rscript

######################################
## ratio WT/Slow from Maslon and m6A # 20220203
######################################

# -------
# Paths |
# -------
rates_Maslon_WT_path <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_data/Elongation_rates_MaslonWT.txt"
rates_Maslon_Slow_path <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_data/Elongation_rates_MaslonSlowSlow.txt"
m6A_WT_path <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_output/3_Normalized_data/meRIP_Ensembl_Normalized_WT_2Kb.txt"
rates_out <- "/media/cc/A/Alicia/Maslon2019_3/Maslon3_output/4_meRIP_Maslon/"

# -----------
# Open data |
# -----------
rates_Maslon_WT <- read.table(rates_Maslon_WT_path,h=T,sep="\t",stringsAsFactors=FALSE)
rates_Maslon_Slow <- read.table(rates_Maslon_Slow_path,h=T,sep="\t",stringsAsFactors=FALSE)
head(rates_Maslon_Slow)
nrow(rates_Maslon_Slow) #1761
nrow(rates_Maslon_WT) #816

m6A_WT <- read.table(m6A_WT_path,h=T,sep="\t",stringsAsFactors=FALSE)
head(m6A_WT)
nrow(m6A_WT) #119007

# ---------------------------------------------------
# Merge common transcripts between WT and Slow/Slow |
# ---------------------------------------------------
rates_Maslon_WT_Slow <- Reduce(function(x,y) merge(x = x, y = y, by = "enst_id", sort = F),
                   list(rates_Maslon_WT[,c("gene_short_name","elongation_b_min","enst_id")],
                        rates_Maslon_Slow[,c("elongation_b_min","enst_id")]))

head(rates_Maslon_WT_Slow)
nrow(rates_Maslon_WT_Slow) #123 !!!!

# ------------otros intentos haciendo merge con gene_short_name: no es fiable

rates_Maslon_WT_Slow2 <- Reduce(function(x,y) merge(x = x, y = y, by = "gene_short_name", sort = F),
                               list(rates_Maslon_WT[,c("gene_short_name","elongation_b_min","enst_id")],
                                    rates_Maslon_Slow[,c("elongation_b_min","enst_id","gene_short_name")]))
head(rates_Maslon_WT_Slow2)
nrow(rates_Maslon_WT_Slow2) #339

rates_Maslon_WT_Slow2[duplicated(rates_Maslon_WT_Slow2$gene_short_name),]
rates_Maslon_WT_Slow$gene_short_name.y == rates_Maslon_WT_Slow$gene_short_name.x

aa <- rates_Maslon_WT[duplicated(rates_Maslon_WT$gene_short_name),]
bb <- rates_Maslon_Slow[duplicated(rates_Maslon_Slow$gene_short_name),]
head(aa)
bb
cc <- Reduce(function(x,y) merge(x = x, y = y, by = "gene_short_name", sort = F),
             list(aa[,c("gene_short_name","elongation_b_min","enst_id")],
                  bb[,c("gene_short_name","elongation_b_min","enst_id")]))

nrow(cc)
rates_Maslon_WT[which(rates_Maslon_WT$gene_short_name == "Slc12a8"), ]
rates_Maslon_Slow[which(rates_Maslon_Slow$gene_short_name == "Slc12a8"), ]
rates_Maslon_WT_Slow2[which(rates_Maslon_WT_Slow2$gene_short_name == "Slc12a8"), ]

# ----------------
# WT/Slow ratios |
# ----------------
rates_Maslon_WT_Slow$WT_Slow_rate = rates_Maslon_WT_Slow$elongation_b_min.x/rates_Maslon_WT_Slow$elongation_b_min.y
head(rates_Maslon_WT_Slow)

# --------------------------------
# Merge m6A WT and Maslon ratios |
# --------------------------------
rates_Maslon_WT_Slow_m6A_WT<- Reduce(function(x,y) merge(x = x, y = y, by.x = "Transcript.stable.ID",by.y = "enst_id", sort = F),
                             list(m6A_WT[,c("Chr", "Start", "End", "Gene.stable.ID","Strand", "Gene_name", "Gene_type",
                                              "meRIP_WT1","meRIP_WT2","meRIP_WT","Transcript.stable.ID")],
                                  rates_Maslon_WT_Slow[,c("elongation_b_min.x","elongation_b_min.y","WT_Slow_rate",
                                                          "gene_short_name","enst_id")]))
head(rates_Maslon_WT_Slow_m6A_WT)
nrow(rates_Maslon_WT_Slow_m6A_WT) #122

# ---------------------
# Remove m6A outliers |
# ---------------------
Q <- quantile(rates_Maslon_WT_Slow_m6A_WT$meRIP_WT, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(rates_Maslon_WT_Slow_m6A_WT$meRIP_WT)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
saved_values_WT <- subset(rates_Maslon_WT_Slow_m6A_WT, rates_Maslon_WT_Slow_m6A_WT$meRIP_WT > low & rates_Maslon_WT_Slow_m6A_WT$meRIP_WT < up)

head(saved_values_WT)
nrow(saved_values_WT) #112

# ------------------
# Correlation test |
# ------------------
cor.test(saved_values_WT$Reads_norm_WT, saved_values_WT$elongation_b_min, method=c("pearson", "kendall", "spearman"))
cor.test(merged$TTseq_WT, merged$RNApol_WT_GB, method=c("pearson"))

# Spearman resulted to be more significant
png(file = paste0(rates_out, "5_Plots/Spearman_MaslonElongationRatio_meRIP_WT_2Kb.png"))
ggscatter(rates_Maslon_WT_Slow_m6A_WT, y = "meRIP_WT", x = "WT_Slow_rate",
          color="lightblue3",shape = 20,
          size=1.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "m6A levels", xlab = "WT/Slow elongation rate ratio Maslon data") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation m6A levels vs. WT/Slow elongation rate ratio \n (Maslon data) - WT (122 genes + 2Kb)")
dev.off()

# Spearman resulted to be more significant
png(file = paste0(rates_out, "5_Plots/Spearman_MaslonElongationRatio_meRIP_WT_2_2Kb.png"))
ggscatter(saved_values_WT, y = "meRIP_WT", x = "WT_Slow_rate",
          color="lightblue3",shape = 20,
          size=1.5,
          add.params = list(color = "lightseagreen", fill = "lightgray"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "m6A levels", xlab = "WT/Slow elongation rate ratio Maslon data") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Spearman correlation m6A levels vs. WT/Slow elongation rate ratio \n (Maslon data) - WT (112 genes + 2Kb)")
dev.off()


# -----------
# Save data |
# -----------
saved_values_WT_ordered <- saved_values_WT[order(saved_values_WT$WT_Slow_rate,decreasing = F),] 
saved_values_WT_ordered <- subset(saved_values_WT_ordered, select=c(Chr, Start, End, gene_short_name, Transcript.stable.ID,
                                                                    Gene.stable.ID, Gene_type,elongation_b_min.x,elongation_b_min.y, 
                                                                    WT_Slow_rate, meRIP_WT))
colnames(saved_values_WT_ordered) <- c("Chr","Start","End","Gene_name","Transcript_ID",
                                       "Gene_ID","Gene_type","Elongation_b_min_WT","Elongation_b_min_SlowSlow",
                                       "WT_SlowSlow_rate", "meRIP_WT")

head(saved_values_WT_ordered)

write.table(saved_values_WT_ordered, 
            file = paste0(rates_out,"meRIP_Normalized_Maslon_WT_noOutliers_2Kb_ordered.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 

rates_Maslon_WT_Slow_m6A_WT

rates_Maslon_WT_Slow_m6A_WT_ordered <- rates_Maslon_WT_Slow_m6A_WT[order(rates_Maslon_WT_Slow_m6A_WT$WT_Slow_rate,decreasing = F),] 
rates_Maslon_WT_Slow_m6A_WT_ordered <- subset(rates_Maslon_WT_Slow_m6A_WT_ordered, select=c(Chr, Start, End, gene_short_name, Transcript.stable.ID,
                                                                    Gene.stable.ID, Gene_type,elongation_b_min.x,elongation_b_min.y, 
                                                                    WT_Slow_rate, meRIP_WT))
colnames(rates_Maslon_WT_Slow_m6A_WT_ordered) <- c("Chr","Start","End","Gene_name","Transcript_ID",
                                       "Gene_ID","Gene_type","Elongation_b_min_WT","Elongation_b_min_SlowSlow",
                                       "WT_SlowSlow_rate", "meRIP_WT")

nrow(rates_Maslon_WT_Slow_m6A_WT_ordered)

write.table(rates_Maslon_WT_Slow_m6A_WT_ordered, 
            file = paste0(rates_out,"meRIP_Normalized_Maslon_WT_withOutliers_2Kb_ordered.txt"),
            quote = F, sep="\t", col.names = T, row.names = F) 


# -----------
# Bar chart | Acabé haciéndolo en excel
# -----------

saved_values_WT2 <- subset(saved_values_WT, select=c(gene_short_name,WT_Slow_rate,meRIP_WT))
head(saved_values_WT2)
saved_values_WT2 <- saved_values_WT2[order(saved_values_WT2$WT_Slow_rate,decreasing = F),] 

library("ggplot2")
library("ggpubr")
library("viridis")
library("RColorBrewer")
library("reshape2")

#saved_values_WT2_melt <- melt(saved_values_WT2, id.vars = "gene_short_name")

ordered<-saved_values_WT2
ordered<-cbind(ordered,tag=1:nrow(saved_values_WT2))
head(ordered)
saved_values_WT2_melt <- melt(ordered, id.vars = c("gene_short_name","tag"))
head(saved_values_WT2_melt)

saved_values_WT2_melt <- saved_values_WT2_melt[order(saved_values_WT2_melt$tag),]
saved_values_WT2_melt$tag<-NULL

df2<-melt(df,id.vars = "tag")
df2<-df2[order(df2$tag),]
df2$tag<-NULL#drop it like it's hot
head(df2)


saved_values_WT2_melt <- saved_values_WT2_melt[order(saved_values_WT2_melt$WT_Slow_rate,decreasing = F),] 
head(saved_values_WT2)
head(aa)
head(saved_values_WT2[,-3])

ggplot(saved_values_WT2_melt, aes(fill = variable, y=value, x=gene_short_name)) + 
  geom_bar(position='stack', stat='identity')


barplot(saved_values_WT4, 
        horiz=TRUE,
        main = "total revenue", 
        names.arg = saved_values_WT3$Gene_name, 
        xlab = "ratioWT/slow", 
        ylab = "Gene", 
        col = c("green","orange"))



# -------------------
# Cumulative curves |
# -------------------

saved_values_WT3 <- saved_values_WT[order(saved_values_WT$WT_Slow_rate,decreasing = F),] 
saved_values_WT4 <- saved_values_WT2[,-1]
saved_values_WT4 <- data.matrix(saved_values_WT4)
nrow(saved_values_WT4)
head(saved_values_WT2)



plot(saved_values_WT3$WT_Slow_rate,
     col="blue",
     main="WT/Slow RNAPII ratio Maslon (112 genes + 2Kb)",
     xlab="WT transcripts",
     ylab="WT/Slow RNAPII ratio",
     mtext("m6A levels", side = 4, line=0),
     lwd = 2,
     type = "l")
lines(saved_values_WT3$meRIP_WT,
      col="grey",
      lwd = 1,
      type = "l")
legend("topright", inset=.02, legend=c("WT/Slow RNAPII ratio", "m6A levels"),
       col=c("blue", "grey"), lty=1, lwd=2, cex=0.8)

