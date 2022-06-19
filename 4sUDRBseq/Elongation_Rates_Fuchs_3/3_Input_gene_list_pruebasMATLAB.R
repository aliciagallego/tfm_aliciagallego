#!/usr/bin/env Rscript

#######################################
## Pruebas para mejorar Matlab script # Tertiles # 20220428
#######################################

# -------
# Paths |
# -------
TTseq_path <- "/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/1.1_Rate_calculation/Elongation_rate_5min_20220401_20Kb_size_processed.txt"
refseq_path <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/4_Input_genes/Input_genes_20Kb.txt")
output <- ("/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/")

# -----------
# Open data |
# -----------
refseq <- read.table(refseq_path,h=T,sep="\t",stringsAsFactors=FALSE)
TTseq <- read.table(TTseq_path,h=T,sep="\t",stringsAsFactors=F,
                    col.names = c("Gene_name","Chr", "Start","End","Exon_number","Strand","TTseq_WT1","TTseq_WT2","TTseq_TKO1","TTseq_TKO2","Size","TTseq_WTmean","TTseq_TKOmean"))

# ---------------------------
# Identifying biased values |
# ---------------------------
TTseq_orderedWT1 <- TTseq[order(TTseq$TTseq_WT1, decreasing = T),]
head(TTseq_orderedWT1, 10L)

TTseq_35 <- TTseq[which(TTseq$TTseq_WT1 == 3.5),]
head(TTseq_35, 10L)

TTseq_05 <- TTseq[which(TTseq$TTseq_WT1 == 0.5 & TTseq$TTseq_WT2 == 0.5),]
head(TTseq_05, 30L)
TTseq_05_2 <- TTseq[which(TTseq$TTseq_TKO1 == 0.5 & TTseq$TTseq_TKO2 == 0.5),]
head(TTseq_05_2, 10L)

# ----------------------------------------------------
# Filter refseq list based to take some genes == 0.5 |
# ----------------------------------------------------
head(refseq)
refseq3 <- refseq[which(refseq$Name == "Mcmdc2" | 
                         refseq$Name == "Gm973" |
                         refseq$Name == "Gpr1" |
                         refseq$Name == "Dner" |
                         refseq$Name == "Sag" |
                         refseq$Name == "B3gat2" |
                          refseq$Name == "Ica1l" |
                          refseq$Name == "Cntnap5b" |
                          refseq$Name == "Pik3c2b" |
                          refseq$Name == "Syt2" |
                          refseq$Name == "Wapl" |
                          refseq$Name == "Prrc2b" |
                          refseq$Name == "Fam168b" |
                          refseq$Name == "Ppp1r12a" |
                          refseq$Name == "Hdac9" |
                          refseq$Name == "Pcdhga6"),]
nrow(refseq2)

# ------------------
# Save merged file |
# ------------------
write.table(refseq3, file = paste0(output,"Input_genes_20Kb_prueba_varios.txt"), 
            quote = F, sep="\t", col.names = T, row.names = F)
