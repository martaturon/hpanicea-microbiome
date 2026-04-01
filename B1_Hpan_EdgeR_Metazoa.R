#Load Libraries
library(edgeR)
library(readxl)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(RColorBrewer) 
library(viridisLite)
library(readr)
library(dplyr)
library(writexl)

##Input tables with isoforms
RSEMtable <- read_tsv("Hpan_final.isoforms.counts.matrix_METAZOA.txt")
counts <- as.data.frame(RSEMtable[,-1])
row.names(counts) <- RSEMtable$Contig

RSEMtableTMM <- read_tsv("Hpan_final.isoform.TMM.EXPR.matrix.txt")
countsTMM <- as.data.frame(RSEMtableTMM[,-1])
row.names(countsTMM) <- RSEMtableTMM$Contig
#Filter Metazoa
countsTMM <- countsTMM[row.names(counts),]
RSEMtableTMM <- RSEMtableTMM %>% 
  filter(Contig %in% row.names(countsTMM))

Infotable <- read_excel("Sample_Information_Hpan.xlsx", sheet = "RNA")
Infotable <- as.data.frame(Infotable)
row.names(Infotable) <- Infotable$Sample_ID

counts <- counts[ , order(colnames(counts)) ]
countsTMM <- countsTMM[ , order(colnames(countsTMM)) ]
Infotable <- Infotable[order(row.names(Infotable)),] 

identical(colnames(counts), row.names(Infotable))
identical(colnames(countsTMM), row.names(Infotable))

#Filter Larvae samples (outliers)
remove <- c("Hp_L1", "Hp_L2")

Infotable <-
  Infotable %>% 
  filter(!Sample_ID %in% remove)

countsTMM <- countsTMM[,rownames(Infotable)]
countsTMM <- countsTMM[rowSums(countsTMM) > 0, ]

identical(colnames(countsTMM), row.names(Infotable))

counts <- counts[,rownames(Infotable)]
counts <- counts[rowSums(counts) > 0, ]

identical(colnames(countsTMM), row.names(Infotable))

#Create subsets to analyse
#####March#####
Infotable_March <-
  Infotable %>% 
  filter(Month == "March")

counts_March <- counts[,row.names(Infotable_March)]
counts_March <- counts_March[rowSums(counts_March) > 0, ]

counts_TMM_March <- countsTMM[,row.names(Infotable_March)]
counts_TMM_March <- counts_TMM_March[rowSums(counts_TMM_March) > 0,]

# CHECK
sampleTotal <- colSums(counts_March)
unique(Infotable_March$Binary)

######March_repro vs NR######
dge <- 
  DGEList(counts = counts_March,  
          group = Infotable_March$Binary,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 25383
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTests(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
# lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
# nrow(lrt.sig)
DA_table <- lrt.sig$table

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_MarchRNR <- all_sigs

######Heatmap
counts_DA <- counts_TMM_March[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
#write_xlsx(DA_table,"Table1_DA_March_RvsNR_Metazoa.xlsx")

######Heatmap March_repro vs NR######
identical(colnames(counts_DA), rownames(Infotable_March))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_March$Status_fine)) 


unique(samp_annot$Repro) 
colors_repro <- c("#365C83","#E3C1CB")
names(colors_repro) <- unique(samp_annot$Repro)

annot_col <- list(Repro=colors_repro)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

####WIthout clustering columns, following Repro vs NR
order <- order(Infotable_March$Binary)
target <- rownames(Infotable_March[order,    ])

pheatmap(log10(counts_DA[,target]+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = F) 

#####April#####
Infotable_April <-
  Infotable %>% 
  filter(Month == "April")

counts_April <- counts[,row.names(Infotable_April)]
counts_April <- counts_April[rowSums(counts_April) > 0, ]

counts_TMM_April <- countsTMM[,row.names(Infotable_April)]
counts_TMM_April <- counts_TMM_April[rowSums(counts_TMM_April) > 0,]

# CHECK
sampleTotal <- colSums(counts_April)
unique(Infotable_April$State)

######April oocytes vs sperm#####
dge <- 
  DGEList(counts = counts_April,  
          group = Infotable_April$State,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 28647 
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTests(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

DA_table <- lrt.sig$table

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_APR_oovssp <- all_sigs
######Heatmap
counts_DA <- counts_TMM_April[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
#write_xlsx(DA_table,"Table2_DA_Apr_RvsNR_Metazoa.xlsx")


######Heatmap April Male vs Fem######
identical(colnames(counts_DA), rownames(Infotable_April))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_April$Status_fine)) 


unique(samp_annot$Repro) 
colors_repro <- c("#AD5A6B","#CDD6AD")
names(colors_repro) <- unique(samp_annot$Repro)

annot_col <- list(Repro=colors_repro)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

####WIthout clustering columns, following Repro vs NR
order <- order(Infotable_April$Sex)
target <- rownames(Infotable_April[order,    ])

pheatmap(log10(counts_DA[,target]+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = F) 

#####May#####
Infotable_May <-
  Infotable %>% 
  filter(Month == "May")

counts_May <- counts[,row.names(Infotable_May)]
counts_May <- counts_May[rowSums(counts_May) > 0, ]

counts_TMM_May <- countsTMM[,row.names(Infotable_May)]
counts_TMM_May <- counts_TMM_May[rowSums(counts_TMM_May) > 0,]

# CHECK
sampleTotal <- colSums(counts_May)
unique(Infotable_May$Sex)

######May_male vs female######
dge <- 
  DGEList(counts = counts_May,  
          group = Infotable_May$Sex,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 26616
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  ## normalmente usariamos ~0 + groups, pero con dos grupos solo no funciona
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTests(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

DA_table <- lrt.sig$table

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_MAY_FvsM <- all_sigs

######Heatmap
counts_DA <- counts_TMM_May[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
#write_xlsx(DA_table,"Table3_DA_May_MALvsFEM_Metazoa.xlsx")


######Heatmap May Male vs Fem######
identical(colnames(counts_DA), rownames(Infotable_May))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_May$Status_fine)) 


unique(samp_annot$Repro) 
colors_repro <- c("#84b8c6","#CDD6AD", "#B46878")
names(colors_repro) <- unique(samp_annot$Repro)

annot_col <- list(Repro=colors_repro)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

####WIthout clustering columns, following Repro vs NR
order <- order(Infotable_May$State)
target <- rownames(Infotable_May[order,    ])

pheatmap(log10(counts_DA[,target]+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = F) 

#####May_without Hp58#####
Infotable_May <-
  Infotable_May %>% 
  filter(!Sample_ID == "Hp_58")

counts_May <- counts[,row.names(Infotable_May)]
counts_May <- counts_May[rowSums(counts_May) > 0, ]

counts_TMM_May <- countsTMM[,row.names(Infotable_May)]
counts_TMM_May <- counts_TMM_May[rowSums(counts_TMM_May) > 0,]

# CHECK
sampleTotal <- colSums(counts_May)
unique(Infotable_May$Sex)

######May no 58 male vs female######
dge <- 
  DGEList(counts = counts_May,  
          group = Infotable_May$Sex,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 26240
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  ## normalmente usariamos ~0 + groups, pero con dos grupos solo no funciona
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTests(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

DA_table <- lrt.sig$table

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_MAY_FvsM <- all_sigs

######Heatmap
counts_DA <- counts_TMM_May[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
#write_xlsx(DA_table,"Table3_DA_May_MALvsFEM_NO58_Metazoa.xlsx")

######Heatmap May Male vs Fem no 58######
identical(colnames(counts_DA), rownames(Infotable_May))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_May$Status_fine)) 


unique(samp_annot$Repro) 
colors_repro <- c("#84b8c6","#CDD6AD", "#B46878")
names(colors_repro) <- unique(samp_annot$Repro)

annot_col <- list(Repro=colors_repro)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

####WIthout clustering columns, following Repro vs NR
order <- order(Infotable_May$State)
target <- rownames(Infotable_May[order,    ])

pheatmap(log10(counts_DA[,target]+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = F) 

#####June#####
Infotable_June <-
  Infotable %>% 
  filter(Month == "June") 

counts_June <- counts[,row.names(Infotable_June)]
counts_June <- counts_June[rowSums(counts_June) > 0, ]

counts_TMM_June <- countsTMM[,row.names(Infotable_June)]
counts_TMM_June <- counts_TMM_June[rowSums(counts_TMM_June) > 0,]

# CHECK
sampleTotal <- colSums(counts_June)
unique(Infotable_June$Status_fine)

######June_RvsNR######
dge <- 
  DGEList(counts = counts_June,  
          group = Infotable_June$Status_fine,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 24908
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  ## normalmente usariamos ~0 + groups, pero con dos grupos solo no funciona
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTests(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

DA_table <- lrt.sig$table

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_JunRNR <- all_sigs

######Heatmap
counts_DA <- counts_TMM_June[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
#write_xlsx(DA_table,"Table4_DA_Jun_RvsNR_Metazoa.xlsx")


######Heatmap June Larve vs NR######
identical(colnames(counts_DA), rownames(Infotable_June))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_June$Status_fine)) 


unique(samp_annot$Repro) 
colors_repro <- c("#365C83","#b494d8")
names(colors_repro) <- unique(samp_annot$Repro)

annot_col <- list(Repro=colors_repro)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

####WIthout clustering columns, following Repro vs NR
order <- order(Infotable_June$Binary)
target <- rownames(Infotable_June[order,    ])

pheatmap(log10(counts_DA[,target]+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = F) 
