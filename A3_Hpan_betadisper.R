# Load Libraries
library(vegan)
library(MASS)
library(ggplot2)
library("ggrepel")
library("gridExtra")
library(dplyr)
library("cowplot")
library(readxl)
library(ggdendro)
library(tidyverse)
library(RColorBrewer)
library(ggnewscale)
library(MoMAColors)

load("Hpan_DataforAnalisis.RData")
load("Colors.RData")

identical(colnames(counts_RA), rownames(Infotable))

##### BETA DIVERSITY ASVs #####
# BETA DIVERSITY 
counts_bc <- vegdist(log2(t(counts_RA)+1), method="bray")  # BEST OUTCOME

#####Stats
#Check for dispersion
disp <-betadisper(counts_bc,Infotable$Date) #or Date
permutest(disp) 

#Permanova
adonis2(counts_bc ~ Date, data = Infotable)
adonis2(counts_bc ~ Date*Status2, data = Infotable)
adonis2(counts_bc ~ Date*Status3, data = Infotable)

#Dendogram
#hc <- hclust(counts_bc, method = "complete")  # "ave"
hc <- hclust(counts_bc, method = "average") 
plot(hc)
plot(hc, hang = -1)

#Edit ggdendogram
ggdendrogram(hc, size = 2)

# Parametric PCoA (of the Bray Curtis dissimilarity) 
counts_bc_Cmds <- cmdscale(counts_bc, k=4, eig = T)  #  Classical (Metric) Multidimensional Scaling #K=4 if we want to plot axes 3 and 4
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

pcoa_shepard <- Shepard(counts_bc, counts_bc_Cmds$points)
cor(pcoa_shepard$x, pcoa_shepard$y)^2 # R^2
plot(pcoa_shepard)

mod <- betadisper(counts_bc , Infotable$Status2, 
                  type = c("median"), bias.adjust = FALSE,
                  sqrt.dist = FALSE, add = FALSE)
boxplot(mod)

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)
head(Infotable)

ggplot(counts_bc_Cmds_points, 
       aes(x = X1, y = X2, colour= Infotable$Month, shape=Infotable$Sex)) + 
  geom_point(size = 5, alpha=0.8) + 
  scale_color_manual(values = Colors_month)+
  labs(color = "Date", shape= "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  #geom_text_repel(aes(label=Infotable$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("Halichondria panicea") + xlab ("PCoA_1 [27.23%]") + ylab ("PCoA_2 [17.28%]") 

#######SUBSETS#####

#####February####
Infotable_Feb <- 
  Infotable %>%
  filter(Month == "FEB")

counts_RA_feb <- counts_RA[,row.names(Infotable_Feb)]
table(rowSums(counts_RA_feb) > 0) 
counts_RA_feb <- counts_RA_feb[rowSums(counts_RA_feb) > 0,]

identical(colnames(counts_RA_feb), rownames(Infotable_Feb))

counts_bc <- vegdist(log2(t(counts_RA_feb)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T , k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

#Check for dispersion
disp <-betadisper(counts_bc,Infotable_Feb$Status2)
permutest(disp) 
#Permanova
adonis2(counts_bc ~ Status2, data = Infotable_Feb)
adonis2(counts_bc ~ Status3, data = Infotable_Feb)

Feb <-
  ggplot(counts_bc_Cmds_points, 
  aes(x = X1, y = X2, colour= Infotable_Feb$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_Feb$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("February") + xlab ("PCoA_1 [34.6%]") + ylab ("PCoA_2 [9.6%]") 

Feb

#####March####
Infotable_march <- 
  Infotable %>%
  filter(Month == "MAR")

counts_RA_march <- counts_RA[,row.names(Infotable_march)]
table(rowSums(counts_RA_march) > 0) 
counts_RA_march <- counts_RA_march[rowSums(counts_RA_march) > 0,]

identical(colnames(counts_RA_march), rownames(Infotable_march))

counts_bc <- vegdist(log2(t(counts_RA_march)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

#Check for dispersion
disp <-betadisper(counts_bc,Infotable_march$Status2)
permutest(disp) 
#Permanova
adonis2(counts_bc ~ Status2, data = Infotable_march)
adonis2(counts_bc ~ Status3, data = Infotable_march)

March <-
  ggplot(counts_bc_Cmds_points, 
  aes(x = X1, y = X2, colour= Infotable_march$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_march$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("March") + xlab ("PCoA_1 [42.42%]") + ylab ("PCoA_2 [16.6%]") 

March

#####April####
Infotable_April <- 
  Infotable %>%
  filter(Month == "APR")

counts_RA_April <- counts_RA[,row.names(Infotable_April)]
table(rowSums(counts_RA_April) > 0) 
counts_RA_April <- counts_RA_April[rowSums(counts_RA_April) > 0,]

identical(colnames(counts_RA_April), rownames(Infotable_April))

counts_bc <- vegdist(log2(t(counts_RA_April)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

#Check for dispersion
disp <-betadisper(counts_bc,Infotable_April$Status2)
permutest(disp) 
#Permanova
adonis2(counts_bc ~ Status2, data = Infotable_April)
adonis2(counts_bc ~ Status3, data = Infotable_April)

April <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X1, y = X2, colour= Infotable_April$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_April$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("April") + xlab ("PCoA_1 [27.7%]") + ylab ("PCoA_2 [18.8%]") 

April

#####May####
Infotable_May <- 
  Infotable %>%
  filter(Month == "MAY")

counts_RA_May <- counts_RA[,row.names(Infotable_May)]
table(rowSums(counts_RA_May) > 0) 
counts_RA_May <- counts_RA_May[rowSums(counts_RA_May) > 0,]

identical(colnames(counts_RA_May), rownames(Infotable_May))

counts_bc <- vegdist(log2(t(counts_RA_May)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

#Check for dispersion
disp <-betadisper(counts_bc,Infotable_May$Status2)
permutest(disp) 
#Permanova
adonis2(counts_bc ~ Status2, data = Infotable_May)
adonis2(counts_bc ~ Status3, data = Infotable_May)

May <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X1, y = X2, colour= Infotable_May$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_May$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("May") + xlab ("PCoA_1 [28.1%]") + ylab ("PCoA_2 [10.2%]") 

May

#####June####
Infotable_June <- 
  Infotable %>%
  filter(Month == "JUN")

counts_RA_June <- counts_RA[,row.names(Infotable_June)]
table(rowSums(counts_RA_June) > 0) 
counts_RA_June <- counts_RA_June[rowSums(counts_RA_June) > 0,]

identical(colnames(counts_RA_June), rownames(Infotable_June))

counts_bc <- vegdist(log2(t(counts_RA_June)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

#Check for dispersion
disp <-betadisper(counts_bc,Infotable_June$Status2)
permutest(disp) 
#Permanova
adonis2(counts_bc ~ Status2, data = Infotable_June)
adonis2(counts_bc ~ Status3, data = Infotable_June)


June <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X1, y = X2, colour= Infotable_June$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_June$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("June") + xlab ("PCoA_1 [36.1%]") + ylab ("PCoA_2 [14.7%]") 

June

#####July####
Infotable_July <- 
  Infotable %>%
  filter(Month == "JUL")

counts_RA_July <- counts_RA[,row.names(Infotable_July)]
table(rowSums(counts_RA_July) > 0) 
counts_RA_July <- counts_RA_July[rowSums(counts_RA_July) > 0,]

identical(colnames(counts_RA_July), rownames(Infotable_July))

counts_bc <- vegdist(log2(t(counts_RA_July)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

July <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X1, y = X2, colour= Infotable_July$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_July$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("July") + xlab ("PCoA_1 [21.8%]") + ylab ("PCoA_2 [14.5%]") 

July

#########ALL###########

library(patchwork) 

(Feb | March | April) / (May | June | July) + plot_layout(guides = 'collect')

########PCoA AXES X3 X4######
#####February####
Infotable_Feb <- 
  Infotable %>%
  filter(Month == "FEB")

counts_RA_feb <- counts_RA[,row.names(Infotable_Feb)]
table(rowSums(counts_RA_feb) > 0) 
counts_RA_feb <- counts_RA_feb[rowSums(counts_RA_feb) > 0,]

identical(colnames(counts_RA_feb), rownames(Infotable_Feb))

counts_bc <- vegdist(log2(t(counts_RA_feb)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T , k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

Feb <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X3, y = X4, colour= Infotable_Feb$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_Feb$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("February") + xlab ("PCoA_3 [8.5%]") + ylab ("PCoA_4 [8%]") 

Feb

#####March####
Infotable_march <- 
  Infotable %>%
  filter(Month == "MAR")

counts_RA_march <- counts_RA[,row.names(Infotable_march)]
table(rowSums(counts_RA_march) > 0) 
counts_RA_march <- counts_RA_march[rowSums(counts_RA_march) > 0,]

identical(colnames(counts_RA_march), rownames(Infotable_march))

counts_bc <- vegdist(log2(t(counts_RA_march)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

March <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X3, y = X4, colour= Infotable_march$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_march$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("March") + xlab ("PCoA_3 [7.57%]") + ylab ("PCoA_4 [7.57%]") 

March

#####April####
Infotable_April <- 
  Infotable %>%
  filter(Month == "APR")

counts_RA_April <- counts_RA[,row.names(Infotable_April)]
table(rowSums(counts_RA_April) > 0) 
counts_RA_April <- counts_RA_April[rowSums(counts_RA_April) > 0,]

identical(colnames(counts_RA_April), rownames(Infotable_April))

counts_bc <- vegdist(log2(t(counts_RA_April)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

April <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X3, y = X4, colour= Infotable_April$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_April$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("April") + xlab ("PCoA_3 [11.1%]") + ylab ("PCoA_4 [8.3%]") 

April

#####May####
Infotable_May <- 
  Infotable %>%
  filter(Month == "MAY")

counts_RA_May <- counts_RA[,row.names(Infotable_May)]
table(rowSums(counts_RA_May) > 0) 
counts_RA_May <- counts_RA_May[rowSums(counts_RA_May) > 0,]

identical(colnames(counts_RA_May), rownames(Infotable_May))

counts_bc <- vegdist(log2(t(counts_RA_May)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

May <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X3, y = X4, colour= Infotable_May$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_May$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("May") + xlab ("PCoA_3 [8.6%]") + ylab ("PCoA_4 [7.5%]") 

May

#####June####
Infotable_June <- 
  Infotable %>%
  filter(Month == "JUN")

counts_RA_June <- counts_RA[,row.names(Infotable_June)]
table(rowSums(counts_RA_June) > 0) 
counts_RA_June <- counts_RA_June[rowSums(counts_RA_June) > 0,]

identical(colnames(counts_RA_June), rownames(Infotable_June))

counts_bc <- vegdist(log2(t(counts_RA_June)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

June <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X3, y = X4, colour= Infotable_June$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_June$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("June") + xlab ("PCoA_3 [11.02%]") + ylab ("PCoA_4 [8.8%]") 

June

#####July####
Infotable_July <- 
  Infotable %>%
  filter(Month == "JUL")

counts_RA_July <- counts_RA[,row.names(Infotable_July)]
table(rowSums(counts_RA_July) > 0) 
counts_RA_July <- counts_RA_July[rowSums(counts_RA_July) > 0,]

identical(colnames(counts_RA_July), rownames(Infotable_July))

counts_bc <- vegdist(log2(t(counts_RA_July)+1), method="bray")  # BEST OUTCOME

counts_bc_Cmds <- cmdscale(counts_bc, eig = T, k=4)  #  Classical (Metric) Multidimensional Scaling
counts_bc_Cmds$eig

Cmds.eig <- round(counts_bc_Cmds$eig,2)
Cmds.eig / sum(Cmds.eig) * 100 # Percentage variance explained by PCoA axes

counts_bc_Cmds_points <- data.frame(counts_bc_Cmds$points)

July <-
  ggplot(counts_bc_Cmds_points, 
         aes(x = X3, y = X4, colour= Infotable_July$Status_fine)) + #shape=Infotable_SP$LevelDepth_1100m, Infotable_SP$Area
  geom_point(size = 3) + 
  scale_color_manual(values = Colors_state)+
  labs(color = "State") +
  theme_light() +
  theme(aspect.ratio=1) +
  geom_text_repel(aes(label=Infotable_July$Sample_ID), size = 3, max.overlaps = 50) +
  ggtitle("July") + xlab ("PCoA_3 [12.5%]") + ylab ("PCoA_4 [9.3%]") 

July

#########ALL###########

library(patchwork) 

(Feb | March | April) / (May | June | July) + plot_layout(guides = 'collect')

 ##Saved 15x8



