library(edgeR)
library(readxl)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

load("Hpan_DataforAnalisis.RData")
load("Colors.RData")

####ALL####
keep <- rowSums((counts_RA) > 0.1) >= 1
table(keep)

counts_RA_fil <- counts_RA[keep,]
counts_fil <- counts[keep,]

taxOR_filt <- taxOR[rownames(counts_RA_fil) , ]

sampleTotal <- colSums(counts_fil)

#####GLM 2 vs 2 ALL######
dge <- DGEList(counts = counts_fil,  
               group = Infotable$Status2,
               lib.size = sampleTotal)

dim(dge)
dge$counts
colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  

dge <- calcNormFactors(dge, method = "TMM") 
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

lrt.sig = topTags(lrt, n = nrow(lrt))  
head(lrt.sig$table)
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
#lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
#nrow(lrt.sig)

####As exact test
y <- DGEList(counts = counts_fil,  
             group = Infotable$Status2,
             lib.size = sampleTotal)

y <- calcNormFactors(y) 
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
plotBCV(y)

et <- exactTest(y, pair=c("Reproductive", "Non_reprodctive"))

summary(DElrt<- decideTestsDGE(et, p=0.05, adjust="BH")) #Change to 0.05
DElrt

et.sig = topTags(et, n = nrow(et))  ## añande el FDR
head(et.sig$table)
et.sig = et.sig[ et.sig$table$FDR < 0.05, ]
nrow(et.sig)
#et.sig = et.sig[ abs(et.sig$table$logFC) > 2, ]
#nrow(et.sig)

## Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% #Change to lrt.sig or et.sig depending on the test
  unlist() %>% 
  unique()

counts_RA_fil_sig <- counts_RA_fil[all_sigs, ]
total <- rowSums(counts_RA_fil_sig)/ncol(counts_RA_fil_sig)
sum(total)
taxOR_sig <- taxOR [ rownames(counts_RA_fil_sig ), ]

range(rowSums(counts_RA_fil_sig))
identical(row.names(taxOR_sig), row.names(counts_RA_fil_sig))
identical(colnames(counts_RA_fil_sig), rownames(Infotable))

####Annotation for columns
samp_annot <- data.frame(row.names = colnames(counts_RA_fil_sig), 
                         group = factor(Infotable$Status2),
                         Month = factor(Infotable$Month))

annot_col <- list(group=Colors_repro, Month=Colors_month)

pheatmap(log10(counts_RA_fil_sig+1), 
         border_color = "NA", 
         fontsize_row = 5, fontsize_col = 5,
         cluster_cols = F,
         color=viridis(100),
         #cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = T) 

###Bubble plot
counts_RA_fil_sig_tb <-
  as_tibble(counts_RA_fil_sig, rownames ='asv')

taxOR_sig_tb <- 
  as_tibble(taxOR_sig, rownames ='asv')

counts_RA_fil_sig_tb_withTaxa <-   
  counts_RA_fil_sig_tb %>% 
  left_join(taxOR_sig_tb)

head(counts_RA_fil_sig_tb_withTaxa )

plot <-
  counts_RA_fil_sig_tb_withTaxa %>% 
  dplyr::select(-c("OTU_ID", "level2taxa","level2taxa_withOTU", "level3taxa_withOTU","level4taxa", "level4taxa_withOTU", "alllevelstaxa")) %>% 
  pivot_longer(cols=-c(asv, level3taxa, alllevelstaxa_withOTU), names_to = "sample", values_to = "Abund" )

plot <- 
  plot %>%
  data.frame(grp = Infotable$Status3[match(Infotable$Sample_ID, plot$sample)])

plot$grp <- factor(plot$grp, levels = c("Non_reprodctive","sperm", "oocytes","embryos","larvae"))

p <- plot %>%  
  ggplot(aes(x= sample, y= alllevelstaxa_withOTU, colour= level3taxa , size = ifelse(Abund == 0, NA, sqrt(Abund)))) +  ## la primera es   2119_2784_16236 Bacteria;Acidobacteria , pero le pone chloroflexi
  geom_point(alpha=0.8) + labs(size = "sqrt(RA)") +
  scale_size_continuous(range=c(0, 9), breaks = waiver()) +
  scale_color_manual(values=Colors_tax, guide = "none")+
  theme(axis.text.y = element_text(size = 4))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

p+ facet_grid(~grp, scales = "free_x", space= "free_x")+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=1.5, linetype="solid"))+
  theme(strip.text.x = element_text(size=8, color="black"))

#####Separate by month####
###FEB#####
Infotable_Feb <- 
  Infotable %>%
  filter(Month == "FEB")

counts_RA_feb <- counts_RA[,row.names(Infotable_Feb)]
counts_RA_feb <- counts_RA_feb[rowSums(counts_RA_feb) > 0,]
counts_feb <- counts[, row.names(Infotable_Feb)]
counts_feb <- counts_feb[rowSums(counts_feb) > 0,]

keep <- rowSums((counts_RA_feb) > 0.05) >= 1  
table(keep)

counts_RA_fil_feb <- counts_RA_feb[keep,]
counts_fil_feb <- counts_feb[keep,]

taxOR_filt_feb <- taxOR[rownames(counts_RA_fil_feb) , ]

sampleTotal <- colSums(counts_fil_feb)

#####GLM 2 vs 2 FEb######
dge <- DGEList(counts = counts_fil_feb,  
               group = Infotable_Feb$Status2,
               lib.size = sampleTotal)

dim(dge)
dge$counts
colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  

dge <- calcNormFactors(dge, method = "TMM") 
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

lrt.sig = topTags(lrt, n = nrow(lrt))  
head(lrt.sig$table)
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

###MARCH#####
Infotable_march <- 
  Infotable %>%
  filter(Month == "MAR")

##NR, oocytes
Infotable_march <- 
  Infotable_march %>%
  filter(Status3 == "Non_reprodctive" | Status3 == "oocytes" ) 

##NR, sperm
Infotable_march <- 
  Infotable_march %>%
  filter(Status3 == "Non_reprodctive" | Status3 == "sperm" ) 

##oocytes vs sperm
Infotable_march <- 
  Infotable_march %>%
  filter(Status3 == "sperm" | Status3 == "oocytes" ) 


counts_RA_march <- counts_RA[,row.names(Infotable_march)]
counts_RA_march <- counts_RA_march[rowSums(counts_RA_march) > 0,]
counts_march <- counts[, row.names(Infotable_march)]
counts_march <- counts_march[rowSums(counts_march) > 0,]

keep <- rowSums((counts_RA_march) > 0.01) >= 1  
table(keep)

counts_RA_fil_march <- counts_RA_march[keep,]
counts_fil_march <- counts_march[keep,]

taxOR_filt_march <- taxOR[rownames(counts_RA_fil_march) , ]

sampleTotal <- colSums(counts_fil_march)

#####GLM 2 vs 2 March######
dge <- DGEList(counts = counts_fil_march,  
               group = Infotable_march$Status3, #Change to Status2 depending on comp.
               lib.size = sampleTotal)

dim(dge)
dge$counts
colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  

dge <- calcNormFactors(dge, method = "TMM") 
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
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

## Plot

all_sigs <-
  lapply(lrt.sig, row.names) %>% #Change to lrt.sig depending on the test
  unlist() %>% 
  unique()

counts_RA_fil_sig <- counts_RA_fil_march[all_sigs, ]
total <- rowSums(counts_RA_fil_sig)/ncol(counts_RA_fil_sig)
sum(total)
taxOR_sig <- taxOR_filt_march [ rownames(counts_RA_fil_sig ), ]

range(rowSums(counts_RA_fil_sig))
identical(row.names(taxOR_sig), row.names(counts_RA_fil_sig))
identical(colnames(counts_RA_fil_sig), rownames(Infotable_march))

####Annotation for columns
samp_annot <- data.frame(row.names = colnames(counts_RA_fil_sig), 
                         group = factor(Infotable_march$Status2),
                         Month = factor(Infotable_march$Month))

colors_group <- c("#0f4b57","#34c06f")
names(colors_group) <- unique(samp_annot$group)

colors_Month <- colorRampPalette(brewer.pal(8, "Set2"))(6)
names(colors_Month) <- unique(samp_annot$Month)


#Color vector
annot_col <- list(group=colors_group, Month=colors_Month)

#Cannot make heatmap with only 1 ASV
pheatmap(log10(counts_RA_fil_sig+1), 
         border_color = "NA",  #It only works when there are no colnames
         fontsize_row = 5, fontsize_col = 5,
         color=viridis(100),
         #cutree_cols = 2,
         annotation_col = samp_annot,
         #annotation_row = row_annot,
         annotation_colors = annot_col,
         show_colnames = T,
         show_rownames = T) 

###Bubble plot
counts_RA_fil_sig_tb <-
  as_tibble(counts_RA_fil_sig, rownames ='asv')

taxOR_sig_tb <- 
  as_tibble(taxOR_sig, rownames ='asv')

counts_RA_fil_sig_tb_withTaxa <-   
  counts_RA_fil_sig_tb %>% 
  left_join(taxOR_sig_tb)

head(counts_RA_fil_sig_tb_withTaxa )

plot <-
  counts_RA_fil_sig_tb_withTaxa %>% 
  dplyr::select(-c("OTU_ID", "level2taxa","level2taxa_withOTU","level3taxa", "level3taxa_withOTU","level4taxa_withOTU", "alllevelstaxa_withOTU")) %>% 
  pivot_longer(cols=-c(asv, level4taxa, alllevelstaxa), names_to = "sample", values_to = "Abund" )

plot <- 
  plot %>%
  data.frame(grp = Infotable_march$Status3[match(Infotable_march$Sample_ID, plot$sample)])

p <- plot %>%   ## ordenamos el phylum
  ggplot(aes(x= sample, y= asv, colour= level4taxa , size = ifelse(Abund == 0, NA, sqrt(Abund)))) +  ## la primera es   2119_2784_16236 Bacteria;Acidobacteria , pero le pone chloroflexi
  geom_point() + labs(size = "sqrt(RA)") +
  scale_size_continuous(range=c(0, 9), breaks = waiver()) +
  #scale_color_manual(values=Colors, guide = guide_legend(ncol=1))+
  theme(axis.text.y = element_text(size = 4))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

p+ facet_grid(~grp, scales = "free_x", space= "free_x")+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=1.5, linetype="solid"))+
  theme(strip.text.x = element_text(size=8, color="black"))

###APRIL#####
Infotable_April <- 
  Infotable %>%
  filter(Month == "APR")

unique(Infotable_April$Status_fine)
#Cannot run comparisons with NR (only 1 sample)

##Change this depending on the comparison
##oocytes vs sperm
Infotable_April <- 
  Infotable_April %>%
  filter(Status_fine == "sperm" | Status_fine == "late oocytes" ) 

counts_RA_April <- counts_RA[,row.names(Infotable_April)]
counts_RA_April <- counts_RA_April[rowSums(counts_RA_April) > 0,]
counts_April <- counts[, row.names(Infotable_April)]
counts_April <- counts_April[rowSums(counts_April) > 0,]

keep <- rowSums((counts_RA_April) > 0.01) >= 1  
table(keep)

counts_RA_fil_April <- counts_RA_April[keep,]
counts_fil_April <- counts_April[keep,]

taxOR_filt_April <- taxOR[rownames(counts_RA_fil_April) , ]

sampleTotal <- colSums(counts_fil_April)

#####GLM 2 vs 2 April######
dge <- DGEList(counts = counts_fil_April,  
               group = Infotable_April$Status_fine,
               lib.size = sampleTotal)

dim(dge)
dge$counts
colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  

dge <- calcNormFactors(dge, method = "TMM") 
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

lrt.sig = topTags(lrt, n = nrow(lrt))  
head(lrt.sig$table)
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
#lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
#nrow(lrt.sig)

#####PLOT######
## Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

counts_RA_fil_sig <- counts_RA_fil_April[all_sigs, ]
total <- rowSums(counts_RA_fil_sig)/ncol(counts_RA_fil_sig)
sum(total)
taxOR_sig <- taxOR_filt_April [ rownames(counts_RA_fil_sig ), ]

###Bubble plot
counts_RA_fil_sig_tb <-
  as_tibble(counts_RA_fil_sig, rownames ='asv')

taxOR_sig_tb <- 
  as_tibble(taxOR_sig, rownames ='asv')

counts_RA_fil_sig_tb_withTaxa <-   
  counts_RA_fil_sig_tb %>% 
  left_join(taxOR_sig_tb)

head(counts_RA_fil_sig_tb_withTaxa )

plot <-
  counts_RA_fil_sig_tb_withTaxa %>% 
  dplyr::select(-c("OTU_ID", "level2taxa","level2taxa_withOTU", "level3taxa_withOTU","level4taxa", "level4taxa_withOTU", "alllevelstaxa")) %>% 
  pivot_longer(cols=-c(asv, level3taxa, alllevelstaxa_withOTU), names_to = "sample", values_to = "Abund" )

plot <- 
  plot %>%
  data.frame(grp = Infotable_April$Status_fine[match(Infotable_April$Sample_ID, plot$sample)])

p <- plot %>%   ## ordenamos el phylum
  ggplot(aes(x= sample, y= alllevelstaxa_withOTU, colour= level3taxa , size = ifelse(Abund == 0, NA, sqrt(Abund)))) +  ## la primera es   2119_2784_16236 Bacteria;Acidobacteria , pero le pone chloroflexi
  geom_point(alpha=0.8) + labs(size = "sqrt(RA)") +
  scale_size_continuous(range=c(0, 9), breaks = waiver()) +
  scale_color_manual(values=Colors_tax, guide = "none")+
  theme(axis.text.y = element_text(size = 4))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

p+ facet_grid(~grp, scales = "free_x", space= "free_x")+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=1.5, linetype="solid"))+
  theme(strip.text.x = element_text(size=8, color="black"))

###May#####
Infotable_May <- 
  Infotable %>%
  filter(Month == "MAY")

unique(Infotable_May$Status_fine)

##Change this depending on the comparison
##late oocytes vs mid embryos
Infotable_May <- 
  Infotable_May %>%
  filter(Status_fine == "late oocytes" | Status_fine == "late oocytes/early embryos/mid embryos" ) 

##late oocytes vs early embryos
Infotable_May <- 
  Infotable_May %>%
  filter(Status_fine == "late oocytes" | Status_fine == "early embryos" ) 

##late oocytes vs sperm
Infotable_May <- 
  Infotable_May %>%
  filter(Status_fine == "late oocytes" | Status_fine == "sperm" ) 

##mid embryos vs early embryos
Infotable_May <- 
  Infotable_May %>%
  filter(Status_fine == "late oocytes/early embryos/mid embryos" | Status_fine == "early embryos" ) 

##mid embryos vs sperm
Infotable_May <- 
  Infotable_May %>%
  filter(Status_fine == "late oocytes/early embryos/mid embryos" | Status_fine == "sperm" )

##early embryos vs sperm
Infotable_May <- 
  Infotable_May %>%
  filter(Status_fine == "early embryos" | Status_fine == "sperm" )

counts_RA_May <- counts_RA[,row.names(Infotable_May)]
counts_RA_May <- counts_RA_May[rowSums(counts_RA_May) > 0,]
counts_May <- counts[, row.names(Infotable_May)]
counts_May <- counts_May[rowSums(counts_May) > 0,]

keep <- rowSums((counts_RA_May) > 0.01) >= 1  
table(keep)

counts_RA_fil_May <- counts_RA_May[keep,]
counts_fil_May <- counts_May[keep,]

taxOR_filt_May <- taxOR[rownames(counts_RA_fil_May) , ]

sampleTotal <- colSums(counts_fil_May)

#####GLM 2 vs 2 May######
dge <- DGEList(counts = counts_fil_May,  
               group = Infotable_May$Status_fine,
               lib.size = sampleTotal)

dim(dge)
dge$counts
colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples)  

dge <- calcNormFactors(dge, method = "TMM") 
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

lrt.sig = topTags(lrt, n = nrow(lrt))  
head(lrt.sig$table)
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
#lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)

## Plot
all_sigs <-
  lapply(lrt.sig, row.names) %>% 
  unlist() %>% 
  unique()

counts_RA_fil_sig <- counts_RA_fil_May[all_sigs, ]
total <- rowSums(counts_RA_fil_sig)/ncol(counts_RA_fil_sig)
sum(total)
taxOR_sig <- taxOR_filt_May [ rownames(counts_RA_fil_sig ), ]

###Bubble plot
counts_RA_fil_sig_tb <-
  as_tibble(counts_RA_fil_sig, rownames ='asv')

taxOR_sig_tb <- 
  as_tibble(taxOR_sig, rownames ='asv')

counts_RA_fil_sig_tb_withTaxa <-   
  counts_RA_fil_sig_tb %>% 
  left_join(taxOR_sig_tb)

head(counts_RA_fil_sig_tb_withTaxa )

plot <-
  counts_RA_fil_sig_tb_withTaxa %>% 
  dplyr::select(-c("OTU_ID", "level2taxa","level2taxa_withOTU", "level3taxa_withOTU","level4taxa", "level4taxa_withOTU", "alllevelstaxa")) %>% 
  pivot_longer(cols=-c(asv, level3taxa, alllevelstaxa_withOTU), names_to = "sample", values_to = "Abund" )

plot <- 
  plot %>%
  data.frame(grp = Infotable_May$Status_fine[match(Infotable_May$Sample_ID, plot$sample)])

p <- plot %>%   ## ordenamos el phylum
  ggplot(aes(x= sample, y= alllevelstaxa_withOTU, colour= level3taxa , size = ifelse(Abund == 0, NA, sqrt(Abund)))) +  ## la primera es   2119_2784_16236 Bacteria;Acidobacteria , pero le pone chloroflexi
  geom_point(alpha=0.8) + labs(size = "sqrt(RA)") +
  scale_size_continuous(range=c(0, 9), breaks = waiver()) +
  scale_color_manual(values=Colors_tax, guide = "none")+
  theme(axis.text.y = element_text(size = 4))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

p+ facet_grid(~grp, scales = "free_x", space= "free_x")+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=1.5, linetype="solid"))+
  theme(strip.text.x = element_text(size=8, color="black"))


###June#####
Infotable_June <- 
  Infotable %>%
  filter(Month == "JUN")

counts_RA_June <- counts_RA[,row.names(Infotable_June)]
counts_RA_June <- counts_RA_June[rowSums(counts_RA_June) > 0,]
counts_June <- counts[, row.names(Infotable_June)]
counts_June <- counts_June[rowSums(counts_June) > 0,]

keep <- rowSums((counts_RA_June) > 0.01) >= 1  
table(keep)

counts_RA_fil_June <- counts_RA_June[keep,]
counts_fil_June <- counts_June[keep,]

taxOR_filt_June <- taxOR[rownames(counts_RA_fil_June) , ]

sampleTotal <- colSums(counts_fil_June)

#####GLM 2 vs 2 June######
dge <- DGEList(counts = counts_fil_June,  
               group = Infotable_June$Status2,
               lib.size = sampleTotal)

dim(dge)
dge$counts
colnames(dge)
plotMDS(dge)

designLM <- model.matrix(~ group, data = dge$samples) 

dge <- calcNormFactors(dge, method = "TMM") 
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

lrt.sig = topTags(lrt, n = nrow(lrt))  
head(lrt.sig$table)
lrt.sig = lrt.sig[ lrt.sig$table$FDR < 0.05, ]
nrow(lrt.sig)
lrt.sig = lrt.sig[ abs(lrt.sig$table$logFC) > 2, ]
nrow(lrt.sig)


