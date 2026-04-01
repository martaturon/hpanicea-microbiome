# Load Libraries
library(tidyverse)
library(readxl) 
library(vegan)

Raw_table <- read_tsv("zotutab_raw.txt")

counts <- as.data.frame(Raw_table)
rownames(counts) <- counts$OTU_ID
counts <- counts[ , -1]

taxOR <- read_excel("Tax_sina.xlsx")
taxOR <- as.data.frame(taxOR)

row.names(taxOR) <- taxOR$OTU_ID

###Remove taxa that belong to "Eukarya" or "Unclassified"
dim(taxOR)
taxOR <- taxOR[!grepl("*Unclassified*", taxOR$L1),]
dim(taxOR) #3478

####Remove taxa that belong to "Chloroplast" or Mithochondria
taxOR <- taxOR[!grepl("*Chloroplast*", taxOR$L4),]
dim(taxOR) #3309

##Apply these filters to our dataset###
dim(counts)
counts <- counts[rownames(taxOR),] # no Eukarya no unclassified
dim(counts) # 3478

identical(row.names(counts), row.names(taxOR))

Infotable_tb <- read_excel("Sample_Information_Hpan.xlsx", sheet = "AMPLICON")


Infotable_tb <-
  Infotable_tb %>%
  filter(Type != "blank") %>%
  filter(Status2 != "MAL")

Infotable <- as.data.frame(Infotable_tb)
rownames(Infotable) <- Infotable$Sample_ID
head(Infotable)

counts <- counts[,row.names(Infotable)]
table(rowSums(counts) > 0) 
counts <- counts[rowSums(counts) > 0,]

identical(sort(colnames(counts)), sort(row.names(Infotable))) #They need to be in the same order to be identical

###Order cols and rows of all the objects
counts <- counts[,order(colnames(counts))]

Infotable <- Infotable[order(row.names(Infotable)), ]

identical(colnames(counts), row.names(Infotable))

taxOR <- taxOR[row.names(counts),] #Filter Tax again since we have removed some samples

# Add taxonomy columns at different levels in Taxonomy file
splitTaxNoProb <- taxOR
splitTaxNoProb <- as.data.frame(splitTaxNoProb)
head(splitTaxNoProb)
dim(splitTaxNoProb)

level2taxa = paste(splitTaxNoProb[,2], splitTaxNoProb[,3], sep = ';') 
head (level2taxa)
taxOR$level2taxa <- level2taxa
head(taxOR)

level2taxa_withOTU <- paste0(level2taxa,  sep="_",splitTaxNoProb[,1])
head(level2taxa_withOTU)
taxOR$level2taxa_withOTU <- level2taxa_withOTU
head (taxOR)

level3taxa = paste(splitTaxNoProb[,2], splitTaxNoProb[,3], splitTaxNoProb[,4], sep = ';') 
head (level3taxa)
taxOR$level3taxa <- level3taxa
head(taxOR)

level3taxa_withOTU <- paste0(level3taxa, sep="_", splitTaxNoProb[,1])
head(level3taxa_withOTU)
taxOR$level3taxa_withOTU <- level3taxa_withOTU
head (taxOR)

level4taxa = paste(splitTaxNoProb[,2], splitTaxNoProb[,3], splitTaxNoProb[,4],  splitTaxNoProb[,5], sep = ';') 
head (level4taxa)
taxOR$level4taxa <- level4taxa
head(taxOR)

level4taxa_withOTU <- paste0(level4taxa, sep="_", splitTaxNoProb[,1])
head(level4taxa_withOTU)
taxOR$level4taxa_withOTU <- level4taxa_withOTU
head (taxOR)

alllevelstaxa =paste(splitTaxNoProb[,2], splitTaxNoProb[,3], splitTaxNoProb[,4], splitTaxNoProb[,5], splitTaxNoProb[,6], splitTaxNoProb[,7], sep = ';') 
head (alllevelstaxa)
taxOR$alllevelstaxa <- alllevelstaxa
head(taxOR)

alllevelstaxa_withOTU <- paste0(alllevelstaxa, sep="_", splitTaxNoProb[,1])
head(alllevelstaxa_withOTU)
taxOR$alllevelstaxa_withOTU <- alllevelstaxa_withOTU
head (taxOR)

taxOR <- taxOR[, -c(2:7)]

### Relative Abundance (RA) and Presence/Absence (PA)
sampleTotals <- colSums(counts)
counts_RA <- t(t(counts) / sampleTotals * 100)
counts_RA <- data.frame(counts_RA)

counts_PA = as.data.frame(ifelse(counts>0,1,0))   


save( counts, Infotable,taxOR,
      counts_RA, counts_PA, 
      file = "Hpan_DataforAnalisis.RData")

#Checking read size
sampleTotals <- colSums(counts)
sampleTotals_PA <- colSums(counts_PA)

myorden1 <- order(colSums(counts), decreasing =T)
sampleTotals[myorden1]
range (sampleTotals)
as.data.frame(sampleTotals)

as.data.frame(sampleTotals_PA)
range (sampleTotals_PA)

## plot the counts
par(mar=c(6, 4, 2, 4))
barplot(sampleTotals[myorden1], las = 2, cex.axis = 0.8, cex.names=0.8, ylab = 'counts',
        ylim=c(0, 60000))

#add the number of sequences
red <- rgb(1, 0, 0, alpha=0.3)
par(new=TRUE)
barplot(sampleTotals_PA[myorden1], las = 2, cex.axis = 0.6, cex.names=0.8, ylab = 'counts',
        ylim=c(0, 30000), col =red, axes=FALSE, alpha = 0.5)

mtext("ASVs",side=4,col="red",line=4) 
axis(4, ylim=c(0,3000), col="red",col.axis="red",las=1)

table( rowSums(counts) > 0)


###Rarefy

min(colSums(counts))

counts_rare <- rrarefy(t(counts), min(colSums(counts))) ## to the minimun value of yoyr colSums(counts)
counts_rare <- as.data.frame(t(counts_rare))
colnames(counts_rare)
colSums(counts_rare)

table(rowSums(counts_rare) > 0)  # there are few ASV that became 0 now, we need to remove them

counts_rare <- counts_rare[ rowSums(counts_rare) > 0, ]
counts_rare <- as.data.frame(counts_rare)
head (counts_rare) 

taxa_rare <- taxOR[ row.names(counts_rare), ] # and tehrefore remove those from the taxa file

sampleTotals_rare <- colSums(counts_rare)

counts_rare_RA <- t(t (counts_rare) / sampleTotals_rare *100)
counts_rare_RA <- as.data.frame(counts_rare_RA)
head(counts_rare)
head(counts_rare_RA)

save(counts_rare, counts_rare_RA, Infotable, taxa_rare, sampleTotals_rare, file = 'Hpan_DataforAnalisis_rare.RData')

##Select COLORS
library(MoMAColors)

#####Colors for L3####
# level 3 taxonomy ----
level3.Counts.RelAb = apply(counts_RA, 2, function(x) tapply(x, taxOR$level3taxa, sum)) 
rowSums(level3.Counts.RelAb)
colSums(counts_RA)

level3.Counts.RelAb.imp  <- (level3.Counts.RelAb)
table(rowSums(level3.Counts.RelAb.imp) > 5) ##Abans ho tenia a 1, mirar si canviar o fer 5
level3.Counts.RelAb.imp <- level3.Counts.RelAb.imp[rowSums(level3.Counts.RelAb.imp) > 5, ]

other3 <- level3.Counts.RelAb[rowSums(level3.Counts.RelAb) <= 5, , drop=F]
colSums(other3)
level3.Counts.RelAb.imp.all <- rbind(level3.Counts.RelAb.imp, "Bacteria;Xothers" = colSums(other3)  )

unique(row.names(level3.Counts.RelAb.imp.all))

Colors_tax <- moma.colors("Klein", 25)
names(Colors_tax) <- unique(row.names(level3.Counts.RelAb.imp.all))

###Colors month###
unique(Infotable$Month)

#Colors_month <- moma.colors("Levine1", 6)
Colors_month <- c("#6B3848", "#818053", "#8B3E50","#D5BB6C", "#474C66","#E0D9B2")
#Feb       #March     #April     #May       #June     #July
names(Colors_month) <- unique(Infotable$Month)

###Colors state5###
unique(Infotable$Status3)

moma.colors("Levine2", 5)
Colors_state5 <- c("#B46878","#365C83","#CDD6AD","#3f9e85","#b494d8")
names(Colors_state5) <- unique(Infotable$Status3)

###Colors state
unique(Infotable$Status_fine)

moma.colors("Levine2", 8)

Colors_state <- c("#E3C1CB","#365C83","#C18292","#CDD6AD","#B46878","#84b8c6","#3f9e85","#b494d8")
#early ooc   #NR    #Mid ooc   #sperm   #late ooc #mid embry #early e  #larvae
names(Colors_state) <- unique(Infotable$Status_fine)

###Colors binary repro
unique(Infotable$Status2)

Colors_repro <- c("#faa281","#2d4262")
names(Colors_repro) <- unique(Infotable$Status2)

save(Colors_tax, Colors_month, Colors_state5, Colors_state,Colors_repro,
     file = "Colors.RData")










