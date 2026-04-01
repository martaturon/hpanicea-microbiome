# Load Libraries
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggdendro)
library(vegan)
library(tidyverse)
library(readxl)
library(viridisLite)
library(MoMAColors)

load("Hpan_DataforAnalisis.RData")
load("Colors.RData")

# dendrogram
counts_bc <- vegdist((t(counts_RA)), method="bray")  
hc <- hclust(counts_bc, method = "complete")  # "ave"

plot(hc)
plot(hc, hang = -1)

ggdendrogram(hc, rotate = T, size = 2)

# level 3 taxonomy ----
level3.Counts.RelAb = apply(counts_RA, 2, function(x) tapply(x, taxOR$level3taxa, sum)) 
rowSums(level3.Counts.RelAb)
colSums(counts_RA)

level3.Counts.RelAb.imp  <- (level3.Counts.RelAb)
table(rowSums(level3.Counts.RelAb.imp) > 5)
level3.Counts.RelAb.imp <- level3.Counts.RelAb.imp[rowSums(level3.Counts.RelAb.imp) > 5, ]

other3 <- level3.Counts.RelAb[rowSums(level3.Counts.RelAb) <= 5, , drop=F]
colSums(other3)
level3.Counts.RelAb.imp.all <- rbind(level3.Counts.RelAb.imp, "Bacteria;Xothers" = colSums(other3)  )


####Separate Can. Halichondribacter
countsRA_ZOTU1 <- counts_RA[row.names(counts_RA) == "Zotu1",]
countsRA_noZOTU1 <- counts_RA[!row.names(counts_RA) == "Zotu1",]
taxOR_noZOTU1 <- taxOR[!row.names(taxOR) == "Zotu1",]

level3.Counts.noZOTU1 = apply(countsRA_noZOTU1, 2, function(x) tapply(x, taxOR_noZOTU1$level3taxa, sum)) 
rowSums(level3.Counts.noZOTU1)
colSums(countsRA_noZOTU1)

level3.Counts.RelAb.imp  <- (level3.Counts.noZOTU1)
table(rowSums(level3.Counts.RelAb.imp) > 20) #To leave only 15 classes
level3.Counts.RelAb.imp <- level3.Counts.RelAb.imp[rowSums(level3.Counts.RelAb.imp) > 20, ]
other3 <- level3.Counts.noZOTU1[rowSums(level3.Counts.noZOTU1) <= 20, , drop=F]
colSums(other3)

level3.Counts.RelAb.imp.all <- rbind(level3.Counts.RelAb.imp, "Bacteria;Xothers" = colSums(other3))
colSums(level3.Counts.RelAb.imp.all)

level3.Counts.RelAb.imp.all_ZOTU1 <- rbind(level3.Counts.RelAb.imp.all, countsRA_ZOTU1)
colSums(level3.Counts.RelAb.imp.all_ZOTU1) #Make sure it is 100

##Change rownames
rownames(level3.Counts.RelAb.imp.all_ZOTU1)[rownames(level3.Counts.RelAb.imp.all_ZOTU1) == "Bacteria;NA;NA"] <- "Bacteria;Xunidentified"
rownames(level3.Counts.RelAb.imp.all_ZOTU1)[rownames(level3.Counts.RelAb.imp.all_ZOTU1) == "Zotu1"] <- "Ca. H. symbioticus"

##L3 BARPLOT#####
level3tax_tb <- 
  as_tibble(level3.Counts.RelAb.imp.all_ZOTU1, rownames ='tax')

plotL3 <-
  level3tax_tb %>% 
  pivot_longer(cols=-c(tax), names_to = "sample", values_to = "Abund" )  ## dejamos los dos niveles por si acasao

plot_dateL3 <-
  plotL3 %>% 
  data.frame(Month = Infotable$Month[match(Infotable$Sample_ID, plotL3$sample)]) %>% 
  data.frame(Repro = Infotable$Status_fine[match(Infotable$Sample_ID, plotL3$sample)]) %>%
  unite(sample_repro, Repro,sample,remove=FALSE)

colnames(plot_dateL3)


##Order##
plot_dateL3$Month <- factor(plot_dateL3$Month, levels = c("FEB","MAR","APR","MAY","JUN","JUL"))
unique(row.names(level3.Counts.RelAb.imp.all_ZOTU1))
Colors_tax <- rev(moma.colors("Rattner", 12))
Colors_tax <- c(Colors_tax,"#8B8989","#251e1f", "#FFFACD")

p <- plot_dateL3 %>%
  #mutate(sample = factor(sample, levels = rev(annot_ord3))) %>%
  ggplot(aes(x=sample_repro, y=Abund, fill=tax)) +
  scale_fill_manual(values=Colors_tax, guide = guide_legend(ncol=1))+
  geom_bar(colour="black",stat="identity") +
  xlab("\nSamples") +
  ylab("Rel. Abundance\n") + 
  theme_bw() +
  theme(axis.text=element_text(size=9),axis.text.x=element_text(angle=90))+
  theme (legend.position="right", 
         legend.key.size = unit(.3,"cm"), #Change the size of the geometric legend elements
         legend.key.width = unit(.3, "cm"))

p

p + facet_grid(~Month, scales = "free_x", space= "free_x")+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=1.5, linetype="solid"))+
  theme(strip.text.x = element_text(size=8, color="black"))

#####With means#####

plot_dateL3_mean <- 
  plot_dateL3 %>%
  group_by(tax, Month, Repro) %>%
  summarise(sum = mean (Abund)) %>%
  unite(Month_repro, Month, Repro, remove=FALSE)

###To put them in the correct order
plot_dateL3_mean$Month <- factor(plot_dateL3_mean$Month, levels = c("FEB","MAR","APR","MAY","JUN","JUL"))

plot_dateL3_mean$Month_repro <- factor(plot_dateL3_mean$Month_repro, 
                                      levels = c("FEB_NR","FEB_early oocytes",
                                                 "MAR_NR","MAR_mid oocytes","MAR_sperm",
                                                 "APR_NR","APR_late oocytes", "APR_sperm",
                                                 "MAY_NR","MAY_late oocytes","MAY_early embryos",
                                                 "MAY_late oocytes/early embryos/mid embryos", "MAY_sperm",
                                                 "JUN_NR","JUN_late embryos/larvae",
                                                 "JUL_NR"))

Colors_tax <- rev(moma.colors("Rattner", 12))
Colors_tax <- c(Colors_tax,"#8B8989","#251e1f", "#FFFACD")

p <- plot_dateL3_mean %>%
  ggplot(aes(x=Month_repro, y=sum, fill=tax)) +
  scale_fill_manual(values=Colors_tax, guide = guide_legend(ncol=1))+
  geom_bar(colour="black",stat="identity") +
  xlab("\nSamples") +
  ylab("Rel. Abundance\n") + 
  theme_bw() +
  theme(axis.text=element_text(size=9),axis.text.x=element_text(angle=90))+
  theme (legend.position="right", 
         legend.key.size = unit(.3,"cm"), #Change the size of the geometric legend elements
         legend.key.width = unit(.3, "cm"))

p + facet_grid(~Month, scales = "free_x", space= "free_x")+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=1.5, linetype="solid"))+
  theme(strip.text.x = element_text(size=8, color="black"))

#Write tables
plot_dateL3_meanX <- 
  plot_dateL3 %>%
  group_by(tax, Month, Repro) %>%
  summarise(sum = mean (Abund)) %>%
  unite(Date_repro, Month, Repro, remove=TRUE)

plot_dateL3_mean_wide <- 
  plot_dateL3_meanX %>%
  pivot_wider(names_from = Date_repro, values_from = sum) 

plot_dateL3_mean_total <- 
  plot_dateL3_mean_wide %>%
  mutate(total = mean(c_across(where(is.numeric)))) %>%
  arrange(desc(total))

#write.csv(plot_dateL3_mean_total , file="Level3_meansZOTU1.csv" )

######LEVEL2
# level 2 taxonomy ----
level2.Counts.RelAb = apply(counts_RA, 2, function(x) tapply(x, taxOR$level2taxa, sum)) 
rowSums(level2.Counts.RelAb)
colSums(counts_RA)

level2.Counts.RelAb.imp  <- (level2.Counts.RelAb)
table(rowSums(level2.Counts.RelAb.imp) > 1)
level2.Counts.RelAb.imp <- level2.Counts.RelAb.imp[rowSums(level2.Counts.RelAb.imp) > 1, ]

other2 <- level2.Counts.RelAb[rowSums(level2.Counts.RelAb) <= 1, , drop=F]
colSums(other2)
level2.Counts.RelAb.imp.all <- rbind(level2.Counts.RelAb.imp, "Bacteria;Xothers" = colSums(other2)  )


##L2 BARPLOT#####
level2tax_tb <- 
  as_tibble(level2.Counts.RelAb.imp.all, rownames ='tax')

plotL2 <-
  level2tax_tb %>% 
  pivot_longer(cols=-c(tax), names_to = "sample", values_to = "Abund" )  ## dejamos los dos niveles por si acasao

plot_dateL2 <-
  plotL2 %>% 
  data.frame(Month = Infotable$Month[match(Infotable$Sample_ID, plotL2$sample)]) %>% 
  data.frame(Repro = Infotable$Status_fine[match(Infotable$Sample_ID, plotL2$sample)])

#####With means#####
plot_dateL2_mean <- 
  plot_dateL2 %>%
  group_by(tax, Month, Repro) %>%
  summarise(sum = mean (Abund)) %>%
  unite(Date_repro, Month, Repro, remove=FALSE)

plot_dateL2_meanX <- 
  plot_dateL2 %>%
  group_by(tax, Month, Repro) %>%
  summarise(sum = mean (Abund)) %>%
  unite(Date_repro, Month, Repro, remove=TRUE)


plot_dateL2_mean_wide <- 
  plot_dateL2_meanX %>%
  pivot_wider(names_from = Date_repro, values_from = sum) 

plot_dateL2_mean_total <- 
  plot_dateL2_mean_wide %>%
  mutate(total = mean(c_across(where(is.numeric)))) %>%
  arrange(desc(total))

#write.csv(plot_dateL2_mean_total , file="Level2_means.txt" )


