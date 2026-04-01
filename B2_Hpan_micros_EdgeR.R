#Load Libraries
library(edgeR)
library(readxl)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(RColorBrewer) 
library(viridisLite)

##Input tables
counts <- read.delim('Hpan_final_isoform_counts_matrix_microsSWPROT.txt')
row.names(counts) <- counts$Feature
counts <- counts[,-1]
countsTMM <- read.delim('Hpan_final_isoform_TMM_matrix_microsSWPROT.txt')
row.names(countsTMM) <- countsTMM$Feature
countsTMM <- countsTMM[,-1]
#counts <- counts_micros
#countsTMM <- tmm_micros

Infotable <- read_excel("Sample_Information_Hpan.xlsx", sheet = "RNA")
Infotable <- as.data.frame(Infotable)
row.names(Infotable) <- Infotable$Sample

#Filter larvae samples (outliers)
remove <- c("Hp_L1", "Hp_L2")

Infotable <-
  Infotable %>% 
  filter(!Sample %in% remove)

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
unique(Infotable_March$Condition)

#####March_repro vs NR######
dge <- 
  DGEList(counts = counts_March,  
          group = Infotable_March$Condition,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 22700 
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  
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

DA_MarchRNR <- all_sigs

######Heatmap
counts_DA <- counts_TMM_March[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table1_DA_March_RvsNR.txt")

identical(colnames(counts_DA), rownames(Infotable_March))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_March$State)) 


unique(samp_annot$Repro) 
colors_repro <- c("#365C83","#AD5A6B")
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

unique(Infotable_April$Condition)
unique(Infotable_April$State)

#####April oocytes vs sperm#####
dge <- 
  DGEList(counts = counts_April,  
          group = Infotable_April$Condition,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 33321 
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  
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

DA_APR_oovssp <- all_sigs
######Heatmap
counts_DA <- counts_TMM_April[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table2_DA_Apr_RvsNR.txt")

identical(colnames(counts_DA), rownames(Infotable_April))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_April$State)) 


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

unique(Infotable_May$Condition2)
unique(Infotable_May$Condition)

#####May_male vs female######
dge <- 
  DGEList(counts = counts_May,  
          group = Infotable_May$Condition,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 2) >= 1 # 22354
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

#fit <- glmFit(dge,designLM)
# Ajustamos el modelo GLM con quasi-likelihood since we have a lot of diversity in the female samples
fit <- glmQLFit(dge, designLM)

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

DA_MAY_FvsM <- all_sigs

######Heatmap
counts_DA <- counts_TMM_May[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table3_DA_May_MALvsFEM_NEW.txt")

identical(colnames(counts_DA), rownames(Infotable_May))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_May$State)) 


unique(samp_annot$Repro) 
colors_repro <- c("#AD5A6B","#CDD6AD")
names(colors_repro) <- unique(samp_annot$Repro)

annot_col <- list(Repro=colors_repro)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         #annotation_col = samp_annot,
         #annotation_row = row_annot,
         #annotation_colors = annot_col,
         show_colnames = T) 

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

unique(Infotable_June$Binary)
unique(Infotable_June$State)

#####June_RvsNR######
dge <- 
  DGEList(counts = counts_June,  
          group = Infotable_June$Binary,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 21845 
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples) 
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
# lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
# nrow(lrt.sig)

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
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table4_DA_Jun_RvsNR.txt")

identical(colnames(counts_DA), rownames(Infotable_June))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Repro = factor(Infotable_June$State)) 


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

######HEATMAP ALL DA REPRO######

DA_ALL <- c(DA_MarchRNR, DA_APR_oovssp, DA_MAY_FvsM, DA_JunRNR)
DA_ALL <- unique(DA_ALL) 

#table
counts_ALL_DA <- countsTMM[DA_ALL,]

counts_ALL_DA <- counts_ALL_DA[,row.names(Infotable)]

identical(colnames(counts_ALL_DA), rownames(Infotable))

samp_annot <- data.frame(row.names = colnames(counts_ALL_DA), 
                         Repro = factor(Infotable$State),
                         Month = factor(Infotable$Month)) 

unique(samp_annot$Repro) 
colors_Repro <- c("#365C83","#AD5A6B","#CDD6AD","#84b8c6","#b494d8")
                   #NR      oocytes    sperm  oocytes/embryos larvae  
names(colors_Repro) <- unique(samp_annot$Repro)

unique(samp_annot$Month) 
colors_Month <- c("#818053", "#6B3848", "#8B3E50", "#D5BB6C")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Repro=colors_Repro, Month= colors_Month)

#Order according to Month and Repro state
counts_ALL_DA_ord <- counts_ALL_DA[ ,c("Hp_16","Hp_17","Hp_26","Hp_20","Hp_27","Hp_30",
                                       "Hp_43","Hp_31","Hp_33","Hp_35","Hp_37","Hp_40","Hp_41","Hp_42",
                                       "Hp_48","Hp_46","Hp_52","Hp_53","Hp_50","Hp_55","Hp_56","Hp_47",
                                       "Hp_61","Hp_66","Hp_69","Hp_68","Hp_71","Hp_74")]

pheatmap(log10(counts_ALL_DA_ord+1), 
         border_color = "NA",  
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         #cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 


#####Males1_AprilMay#####
Infotable_Males1 <-
  Infotable %>% 
  filter(Sex == "MALE") %>% 
  filter(Month == "April" | Month == "May")

counts_Males1 <- counts[,row.names(Infotable_Males1)]
counts_Males1 <- counts_Males1[rowSums(counts_Males1) > 0, ]

counts_TMM_Males1 <- countsTMM[,row.names(Infotable_Males1)]
counts_TMM_Males1 <- counts_TMM_Males1[rowSums(counts_TMM_Males1) > 0,]


# CHECK
sampleTotal <- colSums(counts_Males1)

#GlM
dge <- 
  DGEList(counts = counts_Males1,  
          group = Infotable_Males1$Month,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 5573
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

summary(DElrt<- decideTestsDGE(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
# lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
# nrow(lrt.sig)

DA_table <- lrt.sig$table
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table_DA_males_AprvsMay.txt")

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_Males <- all_sigs

######Heatmap
counts_DA <- counts_TMM_Males1[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table5_DA_Males_AprvsMay.txt")

identical(colnames(counts_DA), rownames(Infotable_Males1))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Month = factor(Infotable_Males1$Month)) 


unique(samp_annot$Month) 
colors_Month <- c("#6B3848", "#8B3E50")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Month=colors_Month)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

#####Females1_MarchApril#####
Infotable_Females1 <-
  Infotable %>% 
  filter(Sex == "FEMALE") %>% 
  filter(Month == "March" | Month == "April")

counts_Females1 <- counts[,row.names(Infotable_Females1)]
counts_Females1 <- counts_Females1[rowSums(counts_Females1) > 0, ]

counts_TMM_Females1 <- countsTMM[,row.names(Infotable_Females1)]
counts_TMM_Females1 <- counts_TMM_Females1[rowSums(counts_TMM_Females1) > 0,]

# CHECK
sampleTotal <- colSums(counts_Females1)

#GLM
dge <- 
  DGEList(counts = counts_Females1,  
          group = Infotable_Females1$Month,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 5662
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTestsDGE(lrt, p=0.05, adjust="BH")) 
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

DA_Fem1 <- all_sigs

######Heatmap
counts_DA <- counts_TMM_Females1[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
write.table(DA_table, row.names = TRUE, sep = "\t", "Table6_DA_Fem_MarchvsApr.txt")

identical(colnames(counts_DA), rownames(Infotable_Females1))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Month = factor(Infotable_Females1$Month)) 


unique(samp_annot$Month) 
colors_Month <- c("#818053","#6B3848")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Month=colors_Month)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

#####Females2_AprilMay#####
Infotable_Females2 <-
  Infotable %>% 
  filter(Sex == "FEMALE") %>% 
  filter(Month == "April" | Month == "May")

counts_Females2 <- counts[,row.names(Infotable_Females2)]
counts_Females2 <- counts_Females2[rowSums(counts_Females2) > 0, ]

counts_TMM_Females2 <- countsTMM[,row.names(Infotable_Females2)]
counts_TMM_Females2 <- counts_TMM_Females2[rowSums(counts_TMM_Females2) > 0,]

# CHECK
sampleTotal <- colSums(counts_Females2)

#GLM
dge <- 
  DGEList(counts = counts_Females2,  
          group = Infotable_Females2$Month,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 5195 
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

summary(DElrt<- decideTestsDGE(lrt, p=0.05, adjust="BH")) 
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

DA_Fem2 <- all_sigs

######Heatmap
counts_DA <- counts_TMM_Females2[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table7_DA_Fem_AprvsMay.txt")

identical(colnames(counts_DA), rownames(Infotable_Females2))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Month = factor(Infotable_Females2$Month)) 


unique(samp_annot$Month) 
colors_Month <- c("#6B3848", "#8B3E50")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Month=colors_Month)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

#####Females3_MayJune#####
Infotable_Females3 <-
  Infotable %>% 
  filter(Sex == "FEMALE") %>% 
  filter(Month == "May" | Month == "June")

counts_Females3 <- counts[,row.names(Infotable_Females3)]
counts_Females3 <- counts_Females3[rowSums(counts_Females3) > 0, ]

counts_TMM_Females3 <- countsTMM[,row.names(Infotable_Females3)]
counts_TMM_Females3 <- counts_TMM_Females3[rowSums(counts_TMM_Females3) > 0,]

# CHECK
sampleTotal <- colSums(counts_Females3)

#GLM
dge <- 
  DGEList(counts = counts_Females3,  
          group = Infotable_Females3$Month,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 4614 
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTestsDGE(lrt, p=0.05, adjust="BH")) 
DElrt

lrt.sig = topTags(lrt, n = nrow(lrt))  ## añande el FDR, LRT = Likelihood Ratio Test
head(lrt.sig$table)
table <- lrt.sig$table
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
# lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
# nrow(lrt.sig)

DA_table <- lrt.sig$table
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table_DA_Females_MaysvsJune.txt")

##Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

DA_Fem3 <- all_sigs

######Heatmap
counts_DA <- counts_TMM_Females3[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table8_DA_Fem_MayvsJun.txt")

identical(colnames(counts_DA), rownames(Infotable_Females3))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Month = factor(Infotable_Females3$Month)) 


unique(samp_annot$Month) 
colors_Month <- c("#8B3E50", "#D5BB6C")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Month=colors_Month)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

#####NR1_MarchJune#####
Infotable_NR1 <-
  Infotable %>% 
  filter(Sex == "NR") %>% 
  filter(Month == "March" | Month == "June")

counts_NR1 <- counts[,row.names(Infotable_NR1)]
counts_NR1 <- counts_NR1[rowSums(counts_NR1) > 0, ]

counts_TMM_NR1 <- countsTMM[,row.names(Infotable_NR1)]
counts_TMM_NR1 <- counts_TMM_NR1[rowSums(counts_TMM_NR1) > 0,]

# CHECK
sampleTotal <- colSums(counts_NR1)

##GLM
dge <- 
  DGEList(counts = counts_NR1,  
          group = Infotable_NR1$Month,
          lib.size = sampleTotal)

dim(dge)
dge$counts
keep <- rowSums(cpm(dge) > 3) >= 2 # 4363 
table(keep)
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # This preforms a TMM normalization (not always necessary)
dge

colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  
#designLM <- model.matrix(~0+group, data=dge$samples) ## con 3 grupos

dge <- estimateGLMCommonDisp(dge,designLM)
dge <- estimateGLMTrendedDisp(dge,designLM)
dge <- estimateGLMTagwiseDisp(dge,designLM)

fit <- glmFit(dge,designLM)

## make contrast here, but we only have one.
lrt <- glmLRT(fit, coef = 2) # compare group2 vs group1
lrt$table
topTags(lrt)

summary(DElrt<- decideTestsDGE(lrt, p=0.05, adjust="BH")) 
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

DA_NR1 <- all_sigs

######Heatmap
counts_DA <- counts_TMM_NR1[all_sigs,]

#Write table
DA_table <- merge(DA_table, counts_DA, by= "row.names")
row.names(DA_table) <- DA_table$Row.names
DA_table <- DA_table[,-1]
#write.table(DA_table, row.names = TRUE, sep = "\t", "Table9_NR_MarvsJun.txt")

identical(colnames(counts_DA), rownames(Infotable_NR1))

samp_annot <- data.frame(row.names = colnames(counts_DA), 
                         Month = factor(Infotable_NR1$Month)) 

unique(samp_annot$Month) 
colors_Month <- c("#818053","#D5BB6C")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Month=colors_Month)

pheatmap(log10(counts_DA+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 

######HEATMAP ALL DA REPRO######

DA_ALL_sex <- c(DA_Males, DA_Fem1, DA_Fem2, DA_Fem3, DA_NR1)
DA_ALL_sex <- unique(DA_ALL_sex) #95

#table
counts_ALL_DA_sex <- countsTMM[DA_ALL_sex,]

counts_ALL_DA_sex <- counts_ALL_DA_sex[,row.names(Infotable)]

identical(colnames(counts_ALL_DA_sex), rownames(Infotable))

samp_annot <- data.frame(row.names = colnames(counts_ALL_DA_sex), 
                         Repro = factor(Infotable$State),
                         Month = factor(Infotable$Month)) 

unique(samp_annot$Repro) 
colors_Repro <- c("#365C83","#AD5A6B","#CDD6AD","#84b8c6","#b494d8")
                    #NR      oocytes    sperm  oocytes/embryos larvae  
names(colors_Repro) <- unique(samp_annot$Repro)

unique(samp_annot$Month) 
colors_Month <- c("#818053", "#6B3848", "#8B3E50", "#D5BB6C")
names(colors_Month) <- unique(samp_annot$Month)

annot_col <- list(Repro=colors_Repro, Month= colors_Month)

#Order according to Month and Repro state
counts_ALL_DA_sex_ord <- counts_ALL_DA_sex[ ,c("Hp_20","Hp_26","Hp_27","Hp_30","Hp_31","Hp_33","Hp_35","Hp_43",
                                               "Hp_48","Hp_46","Hp_52","Hp_53",
                                               "Hp_68","Hp_71","Hp_74",
                                               "Hp_37","Hp_40","Hp_41","Hp_42","Hp_47","Hp_50","Hp_55","Hp_56",
                                               "Hp_16","Hp_17",
                                               "Hp_61","Hp_66","Hp_69")]

pheatmap(log10(counts_ALL_DA_sex_ord+1), 
         border_color = "NA",  
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         cluster_cols = F,
         #cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T) 
