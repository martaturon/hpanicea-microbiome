# Load Libraries
library(vegan)
library(ggplot2)
library(RColorBrewer) 
library(patchwork)
library(colorspace)
library(ggnewscale)
library(tidyr)

load("Hpan_DataforAnalisis.RData")
load("Hpan_DataforAnalisis_rare.RData")
load("Colors.RData")

specNumber_rare <- specnumber(counts_rare, MARGIN = 2)  
Infotable$specNumber_rare <- specNumber_rare  # This is the stimated Richness

shannonH <- diversity(t(counts_rare), index = "shannon")   ## this is diversity
summary(shannonH)
Infotable$shannonH <- shannonH
invSimp <- diversity(t(counts_rare), index = "invsimpson")   ## this is another diversity measure
Infotable$invSimp <- invSimp

Infotable <-
  Infotable %>%
  unite(Month_Stage, Month, Status_fine, remove=F)

ggplot(Infotable, aes(x = Month_Stage, y = shannonH)) +  
  geom_point(size = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Infotable$Month <- factor(Infotable$Month, levels = c("FEB","MAR","APR","MAY","JUN","JUL"))
Infotable$Month_Stage <- factor(Infotable$Month_Stage, 
                               levels = c("FEB_NR","FEB_early oocytes",
                                          "MAR_NR","MAR_mid oocytes","MAR_sperm",
                                          "APR_NR","APR_late oocytes", "APR_sperm",
                                          "MAY_NR","MAY_late oocytes","MAY_late oocytes/early embryos/mid embryos",
                                          "MAY_early embryos", "MAY_sperm",
                                          "JUN_NR","JUN_late embryos/larvae",
                                          "JUL_NR"))

Infotable$Status_fine <- factor(Infotable$Status_fine,
                                levels = c("NR","early oocytes","mid oocytes","late oocytes",
                                           "late oocytes/early embryos/mid embryos","early embryos",
                                           "late embryos/larvae", "sperm"))
####PLOT#####
Hplot <- ggplot(Infotable, aes(x=Month_Stage, y=shannonH )) +
  geom_boxplot(
    aes(fill = Status_fine), 
    size=0.5, alpha=0.8)+
  scale_fill_manual(values=Colors_state)+
  labs(fill = "State") +
  geom_vline(xintercept = c(2.5,5.5,8.5, 13.5, 15.5), linetype = "dashed", color= "grey")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_jitter(size=2, position=position_jitter(width=.25), alpha=0.8, show.legend = NA)
Hplot

#Months mean Shannon
ggplot(Infotable, aes(x=Month, y=shannonH))+
  geom_boxplot(
    aes(color = Month,
        fill= after_scale(desaturate(lighten(color, .6), .6))), 
    size=1.4)+
  scale_color_manual(values=Colors_month)+
  geom_jitter(aes(colour=Month), size=4, position=position_jitter(width=.25), alpha=0.75, show.legend = NA)+
  theme_bw()+
  labs(x="",
       y= "Shannon")+
  theme(legend.position = "none")

Richness <- ggplot(Infotable, aes(x=Month_Stage, y=specNumber_rare )) +
  geom_boxplot(
    aes(fill = Status_fine), 
    size=0.5, alpha=0.8)+
  scale_fill_manual(values=Colors_state)+
  geom_vline(xintercept = c(2.5,5.5,8.5, 13.5, 15.5), linetype = "dashed", color= "grey")+
  labs(x="",
       y= "Richness",
       fill="State")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_jitter(size=2, position=position_jitter(width=.25), alpha=0.8, show.legend = NA)
Richness

#Months mean Richness
ggplot(Infotable, aes(x=Month, y=specNumber_rare))+
  geom_boxplot(
    aes(color = Month,
        fill= after_scale(desaturate(lighten(color, .6), .6))), 
    size=1.4)+
  scale_color_manual(values=Colors_month)+
  geom_jitter(aes(colour=Month), size=4, position=position_jitter(width=.25), alpha=0.75, show.legend = NA)+
  theme_bw()+
  labs(x="",
       y= "Richness")+
  theme(legend.position = "none")

InvSim <- ggplot(Infotable, aes(x=Month_Stage, y=invSimp )) +
  geom_boxplot(
    aes(fill = Status_fine), 
    size=0.5, alpha=0.8)+
  scale_fill_manual(values=Colors_state)+
  labs(x="",
       y= "InvSim",
       fill="State")+
  geom_vline(xintercept = c(2.5,5.5,9.5, 13.5, 15.5), linetype = "dashed", color= "grey")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_jitter(size=2, position=position_jitter(width=.25), alpha=0.8, show.legend = NA)
InvSim

#######Stats######
library(effects)
##first test normality
shapiro.test(Infotable[,"shannonH"]) ###No significatiu, so normality

res.aov_Month <- aov(shannonH ~ Month, data = Infotable)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
tukey.test <- TukeyHSD(res.aov_Month) ##Pair-wise comparisons to see where are the differences
plot(tukey.test)

shannonH_model <- lm(shannonH ~ Month, data = Infotable)
anova(shannonH_model)

TukeyHSD(aov(shannonH ~ Month, data = Infotable))
plot(allEffects(shannonH_model))

shannonH_model <- lm(shannonH ~ Status_fine, data = Infotable)
anova(shannonH_model)

TukeyHSD(aov(shannonH ~ Status_fine, data = Infotable))
plot(allEffects(shannonH_model))

###richness
specNumber_rare_model <- lm(specNumber_rare ~ Month, data = Infotable)
anova(specNumber_rare_model)

TukeyHSD(aov(specNumber_rare ~ Month, data = Infotable))
plot(allEffects(specNumber_rare_model))
###

res.aov_Month <- aov(shannonH ~ Month_Stage, data = Infotable)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
tukey.test <- TukeyHSD(res.aov_Month) ##Pair-wise comparisons to see where are the differences

# Extract significant comparisons
significant_comparisons <- tukey.test$Month_Stage[, "p adj"] < 0.05
significant_pairs <- tukey.test$Month_Stage[significant_comparisons, ]

###Subsets
###FEB#####
Infotable_Feb <- 
  Infotable %>%
  filter(Month == "FEB")

res.aov_Month <- aov(shannonH ~ Month_Stage, data = Infotable_Feb)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
#TukeyHSD(res.aov_Month) 
res.aov_Month <- aov(specNumber_rare ~ Month_Stage, data = Infotable_Feb)
summary(res.aov_Month)

###MARCH#####
Infotable_march <- 
  Infotable %>%
  filter(Month == "MAR")

res.aov_Month <- aov(shannonH ~ Month_Stage, data = Infotable_march)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
#TukeyHSD(res.aov_Month) 
res.aov_Month <- aov(specNumber_rare ~ Month_Stage, data = Infotable_march)
summary(res.aov_Month)

###APRIL#####
Infotable_April <- 
  Infotable %>%
  filter(Month == "APR")

res.aov_Month <- aov(shannonH ~ Month_Stage, data = Infotable_April)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
#TukeyHSD(res.aov_Month) 
res.aov_Month <- aov(specNumber_rare ~ Month_Stage, data = Infotable_April)
summary(res.aov_Month)


###May#####
#We need to take this specific subset only comparing oocytes vs embrio
Infotable_May <- 
  Infotable %>%
  filter(Month == "MAY") #%>%
  #filter(Status_fine == "oocytes" | Status_fine =="sperm") 

res.aov_Month <- aov(specNumber_rare ~ Status_fine, data = Infotable_May)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
stats <- TukeyHSD(res.aov_Month) 

res.aov_Month <- aov(shannonH ~ Status_fine, data = Infotable_May)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
TukeyHSD(res.aov_Month) 

####June###
Infotable_June <- 
  Infotable %>%
  filter(Month == "JUN")

res.aov_Month <- aov(shannonH ~ Month_Stage, data = Infotable_June)
summary(res.aov_Month) ##Significant diff depending on Month for Shannon div
#TukeyHSD(res.aov_Month) 
res.aov_Month <- aov(specNumber_rare ~ Month_Stage, data = Infotable_June)
summary(res.aov_Month)



