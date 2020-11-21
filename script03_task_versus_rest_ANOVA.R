library(car); library(ez); library(tidyverse)

#### OVERALL MODEL -------------------------------------------------------------

DATA2 <- read.csv("./data_fac_anova_25_11182020.csv", 
                  header = TRUE,na.strings=c("","NA","na"),
                  stringsAsFactors = TRUE)

DATA2$subID <- factor(DATA2$subID)

DATA2$time <- fct_recode(DATA2$block, Rest1 = "A", Task = "B", Rest2 = "C")
DATA2$time <- fct_relevel(DATA2$time, "Rest1", "Task", "Rest2")

DATA2$method <- fct_recode(DATA2$method, FOOOF = "fooof", LMER = "lmer")
DATA2$method <- fct_relevel(DATA2$method, "LMER", "FOOOF")


DATA2$channel <- fct_recode(DATA2$channel, Frontal = "frontal", Central = "central",
                           Parietal = "parietal", Occipital = "occipital")
DATA2$channel <- fct_relevel(DATA2$channel, "Frontal", "Central", "Parietal", "Occipital")


# Converting all of the FOOOF exponents to be negative to match the LMER slopes:
DATA2$slope <- ifelse(DATA2$method=="FOOOF", -1*DATA2$slope, DATA2$slope)
head(DATA2)

# Implementing the full ANOVA model
overall_model<-ezANOVA(data=DATA2, dv=.(slope), wid=.(subID), 
                       within = .(method, time, channel), detailed = TRUE, type = 3)

overall_model

t.test(DATA2$slope~DATA2$method)
# The Significant Method x Time x Channel Interaction justifies breaking up the 
# data into smaller ANOVAs. 


# Plot of Slope as a function of Task and Region -------------------------------
head(DATA2)

LMER <- subset(DATA2, method=="LMER")
LMER_model<-ezANOVA(data=LMER, dv=.(slope), wid=.(subID), 
                       within = .(time, channel), detailed = TRUE, type = 3)

LMER_model


ggplot(data = DATA2[DATA2$method=="LMER",], 
       mapping = aes(x = channel, y = slope)) +
  geom_point(aes(fill=time), pch=21, size=2, position=position_dodge(width=0.75))+
  #geom_line(aes(group = subID))+
  geom_boxplot(aes(fill=time), alpha = .6, notch=FALSE, 
               col="black", lwd=1, outlier.shape=NA)+
  scale_x_discrete(name = "Region") +
  scale_y_continuous(name = "Spectral Slope", limits=c(-2.5,0.5)) +
  scale_fill_grey()+
  theme_bw()+
  theme(axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        legend.position = "none")



FOOOF <- subset(DATA2, method=="FOOOF")
FOOOF_model<-ezANOVA(data=FOOOF, dv=.(slope), wid=.(subID), 
                    within = .(time, channel), detailed = TRUE, type = 3)

FOOOF_model


ggplot(data = DATA2[DATA2$method=="FOOOF",], 
       mapping = aes(x = channel, y = slope)) +
  geom_point(aes(fill=time), pch=21, size=2, position=position_dodge(width=0.75))+
  #geom_line(aes(group = subID))+
  geom_boxplot(aes(fill=time), alpha = .6, notch=FALSE, 
               col="black", lwd=1, outlier.shape=NA)+
  scale_x_discrete(name = "Region") +
  scale_y_continuous(name = "Aperiodic Slope", limits=c(-2.5,0.5)) +
  scale_fill_grey()+
  theme_bw()+
  theme(axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        legend.position = "none")




##### Frontal ANOVA Model -------------------------------------------------
FRONTAL <- subset(DATA2, channel =="Frontal")

frontal_model<-ezANOVA(data=FRONTAL, dv=.(slope), wid=.(subID), 
                       within = .(method, time), detailed = TRUE, type = 3)

frontal_model

ggplot(data = FRONTAL, 
       mapping = aes(x = time, y = slope)) +
  geom_point(aes(fill=method), pch=21, size=2)+
  geom_line(aes(group = subID, col=method))+
  geom_boxplot(aes(fill=method, group=time), alpha = .6, notch=FALSE, 
               col="black", lwd=1, outlier.shape=NA)+
  facet_wrap(~method)+
  scale_x_discrete(name = "Experimental Session") +
  scale_y_continuous(name = "PSD Slope") +
  theme_bw()+
  theme(axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        legend.position = "none")


#### Central ANOVA Model --------------------------------------------------
CENTRAL <- subset(DATA2, channel =="Central")

central_model<-ezANOVA(data=CENTRAL, dv=.(slope), wid=.(subID), 
                       within = .(method, time), detailed = TRUE, type = 3)

central_model

#### ANOVA Post hoc testing within lmer method ----------------------------------
CENTRAL_lmer <- subset(CENTRAL, CENTRAL$method == "LMER")
head(CENTRAL_lmer)
t.test(CENTRAL_lmer$slope[CENTRAL_lmer$time=="Task"],
       CENTRAL_lmer$slope[CENTRAL_lmer$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

t.test(CENTRAL_lmer$slope[CENTRAL_lmer$time=="Task"],
       CENTRAL_lmer$slope[CENTRAL_lmer$time=="Rest1"],
       paired=TRUE, var.equal = FALSE)


t.test(CENTRAL_lmer$slope[CENTRAL_lmer$time=="Rest1"],
       CENTRAL_lmer$slope[CENTRAL_lmer$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)


#### ANOVA Post hoc testing within fooof ---------------------------------------

CENTRAL_fooof <- subset(CENTRAL, CENTRAL$method == "FOOOF")

t.test(CENTRAL_fooof$slope[CENTRAL_fooof$time=="Task"],
       CENTRAL_fooof$slope[CENTRAL_fooof$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

t.test(CENTRAL_fooof$slope[CENTRAL_fooof$time=="Task"],
       CENTRAL_fooof$slope[CENTRAL_fooof$time=="Rest1"],
       paired=TRUE, var.equal = FALSE)


t.test(CENTRAL_fooof$slope[CENTRAL_fooof$time=="Rest1"],
       CENTRAL_fooof$slope[CENTRAL_fooof$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)


ggplot(data = CENTRAL, 
       mapping = aes(x = time, y = slope)) +
  geom_point(aes(fill=method), pch=21, size=2)+
  geom_line(aes(group = subID, col=method))+
  geom_boxplot(aes(fill=method, group=time), alpha = .6, notch=FALSE, 
               col="black", lwd=1, outlier.shape=NA)+
  facet_wrap(~method)+
  scale_x_discrete(name = "Experimental Session") +
  scale_y_continuous(name = "PSD Slope") +
  theme_bw()+
  theme(axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        legend.position = "none")


##### Parietal ANOVA Model -------------------------------------------
PARIETAL <- subset(DATA2, channel =="Parietal")

parietal_model<-ezANOVA(data=PARIETAL, dv=.(slope), wid=.(subID), 
                        within = .(method, time), detailed = TRUE, type = 3)

parietal_model


ggplot(data = PARIETAL, 
       mapping = aes(x = time, y = slope)) +
  geom_point(aes(fill=method), pch=21, size=2)+
  geom_line(aes(group = subID, col=method))+
  geom_boxplot(aes(fill=method, group=time), alpha = .6, notch=FALSE, 
               col="black", lwd=1, outlier.shape=NA)+
  facet_wrap(~method)+
  scale_x_discrete(name = "Experimental Session") +
  scale_y_continuous(name = "PSD Slope") +
  theme_bw()+
  theme(axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        legend.position = "none")


### Parietal post hoc testing within lmer -------------------------------
PARIETAL_lmer <- subset(PARIETAL, PARIETAL$method == "LMER")

t.test(PARIETAL_lmer$slope[PARIETAL_lmer$time=="Task"],
       PARIETAL_lmer$slope[PARIETAL_lmer$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

t.test(PARIETAL_lmer$slope[PARIETAL_lmer$time=="Task"],
       PARIETAL_lmer$slope[PARIETAL_lmer$time=="Rest1"],
       paired=TRUE, var.equal = FALSE)


t.test(PARIETAL_lmer$slope[PARIETAL_lmer$time=="Rest1"],
       PARIETAL_lmer$slope[PARIETAL_lmer$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)


#### Parietal Post hoc testing within fooof -----------------------------

PARIETAL_fooof <- subset(PARIETAL, PARIETAL$method == "FOOOF")

t.test(PARIETAL_fooof$slope[PARIETAL_fooof$time=="Task"],
       PARIETAL_fooof$slope[PARIETAL_fooof$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

t.test(PARIETAL_fooof$slope[PARIETAL_fooof$time=="Task"],
       PARIETAL_fooof$slope[PARIETAL_fooof$time=="Rest1"],
       paired=TRUE, var.equal = FALSE)


t.test(PARIETAL_fooof$slope[PARIETAL_fooof$time=="Rest1"],
       PARIETAL_fooof$slope[PARIETAL_fooof$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

#### Occipital ANOVA Model ----------------------------------------------
OCCIPITAL <- subset(DATA2, channel =="Occipital")

occipital_model<-ezANOVA(data=OCCIPITAL, dv=.(slope), wid=.(subID), 
                         within = .(method, time), detailed = TRUE, type = 3)

occipital_model

ggplot(data = OCCIPITAL, 
       mapping = aes(x = time, y = slope)) +
  geom_point(aes(fill=method), pch=21, size=2)+
  geom_line(aes(group = subID, col=method))+
  geom_boxplot(aes(fill=method, group=time), alpha = .6, notch=FALSE, 
               col="black", lwd=1, outlier.shape=NA)+
  facet_wrap(~method)+
  scale_x_discrete(name = "Experimental Session") +
  scale_y_continuous(name = "PSD Slope") +
  theme_bw()+
  theme(axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        legend.position = "none")

### Occipital post hoc testing within lmer -------------------------------
OCCIPITAL_lmer <- subset(OCCIPITAL, method == "LMER")

t.test(OCCIPITAL_lmer$slope[OCCIPITAL_lmer$time=="Task"],
       OCCIPITAL_lmer$slope[OCCIPITAL_lmer$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

t.test(OCCIPITAL_lmer$slope[OCCIPITAL_lmer$time=="Task"],
       OCCIPITAL_lmer$slope[OCCIPITAL_lmer$time=="Rest1"],
       paired=TRUE, var.equal = FALSE)


t.test(OCCIPITAL_lmer$slope[OCCIPITAL_lmer$time=="Rest1"],
       OCCIPITAL_lmer$slope[OCCIPITAL_lmer$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)


#### Parietal Post hoc testing within fooof -----------------------------

OCCIPITAL_fooof <- subset(OCCIPITAL, method == "FOOOF")

t.test(OCCIPITAL_fooof$slope[OCCIPITAL_fooof$time=="Task"],
       OCCIPITAL_fooof$slope[OCCIPITAL_fooof$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)

t.test(OCCIPITAL_fooof$slope[OCCIPITAL_fooof$time=="Task"],
       OCCIPITAL_fooof$slope[OCCIPITAL_fooof$time=="Rest1"],
       paired=TRUE, var.equal = FALSE)


t.test(OCCIPITAL_fooof$slope[OCCIPITAL_fooof$time=="Rest1"],
       OCCIPITAL_fooof$slope[OCCIPITAL_fooof$time=="Rest2"],
       paired=TRUE, var.equal = FALSE)
