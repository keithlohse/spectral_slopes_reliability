##### Opening libraries ----------------------------------------------
library("car"); library("compute.es"); library("tidyverse"); library("multcomp");
library("pastecs"); library("lme4"); library("lmerTest");
library ("ICC"); library("psych"); library("ez"); 




# let's see what is in the data folder
list.files()


#### Log transformation -------------------------------------------------

DATA <- read.csv("./data_EEG_MASTER4.csv", 
                 header = TRUE,na.strings=c("","NA","na"),
                 stringsAsFactors = TRUE)
head(DATA)



DATA$lgFrontal<-log(DATA$Frontal)
DATA$lgCentral<-log(DATA$Central)
DATA$lgParietal<-log(DATA$Parietal)
DATA$lgOccipital<-log(DATA$Occipital)
DATA$lgHz<-log(DATA$Hz)



head(DATA)

## Filter data to frequecy 2-25 hz

FILTERED<-subset(DATA, Hz>=2) # Lose everything below 1
FILTERED<-subset(FILTERED, Hz<=25) # Lose everything above 35

head(FILTERED)
summary(as.factor(FILTERED$Hz))  ## Shows us the number of observations for 
## each frequency in our data


# FIGURE 1A: PSD Showing the Raw Data across positions
DATA_AVE <- FILTERED %>%
            group_by(block, Hz) %>%
            summarize(lgHz = lgHz[1],
                      ave_central = mean(Central, na.rm=TRUE),
                      ave_lg_central = mean(lgCentral, na.rm=TRUE))

head(DATA_AVE)

ggplot(data = DATA_AVE[DATA_AVE$block=="pretest",], 
       mapping = aes(x = Hz, y = ave_central)) +
  scale_x_continuous(name = "Frequency (Hz)") +
  scale_y_continuous(name = "Mean Central Power (uV^2)")+
  annotate("rect", xmin=2.93, xmax=4, ymin=0, ymax=5, 
           alpha=0.2, fill="grey", col="black")+
  annotate("rect", xmin=4, xmax=8, ymin=0, ymax=3.2, 
           alpha=0.2, fill="blue", col="black")+
  annotate("rect", xmin=8, xmax=13, ymin=0, ymax=2.5, 
           alpha=0.2, fill="red", col="black")+
  annotate("rect", xmin=13, xmax=25, ymin=0, ymax=1.9, 
           alpha=0.2, fill="green", col="black")+
  geom_line(col="black", lwd=1) +
  theme_classic()+
  theme(axis.text=element_text(size=16, colour="black"), 
        strip.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom")


# FIGURE 1B: Spectral Slope
ggplot(data = DATA_AVE[DATA_AVE$block=="pretest",], 
       mapping = aes(x = lgHz, y = ave_lg_central)) +
  scale_x_continuous(name = "Log Frequency (ln(Hz))") +
  scale_y_continuous(name = "Log Mean Central Power (ln(uV^2))",
                     limits=c(-2,2))+
  geom_line(col="black", lwd=1) +
  stat_smooth(col="red", method="lm", se=FALSE, lty=2)+
  theme_classic()+
  theme(axis.text=element_text(size=16, colour="black"), 
        strip.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom")

# FIGURE 1C: Aperiodic Slope
ggplot(data = DATA_AVE[DATA_AVE$block=="pretest",], 
       mapping = aes(x = Hz, y = ave_lg_central)) +
  scale_x_continuous(name = "Frequency (Hz)") +
  scale_y_continuous(name = "Log Mean Central Power (ln(uV^2))",
                     limits=c(-2,2))+
  geom_line(col="black", lwd=1) +
  #stat_smooth(col="red", method="lm", se=FALSE, lty=2)+
  theme_classic()+
  theme(axis.text=element_text(size=16, colour="black"), 
        strip.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        legend.position = "bottom")


# Figure 1D: Rest versus Task Activity in the Current Data ---------------------
head(DATA)
summary(DATA$time)
DATA_AVE <- FILTERED %>%
  group_by(time, Hz) %>%
  summarize(lgHz = lgHz[1],
            ave_central = mean(Central, na.rm=TRUE),
            ave_lg_central = mean(lgCentral, na.rm=TRUE))

summary(DATA_AVE$time)
levels(DATA_AVE$time) <- c("On Task", "At Rest")
DATA_AVE$time<-factor(DATA_AVE$time, levels=c("At Rest", "On Task"))

head(DATA_AVE)

ggplot(data = DATA_AVE, 
       mapping = aes(x = Hz, y = ave_central)) +
  scale_x_continuous(name = "Frequency (Hz)") +
  scale_y_continuous(name = "Mean Central Power (uV^2)")+
  geom_line(aes(col=time, lty=time), lwd=1) +
  scale_color_manual(values=c("black", "grey40"))+
  labs(lty = "Condition", col="Condition")+
  theme_classic()+
  theme(axis.text=element_text(size=16, colour="black"), 
        strip.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = c(0.7,0.7),
        legend.background = element_rect(fill="grey90",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))





#### Generating lmer slopes and intercepts ----------------------------------
head(FILTERED)


## generating slopes for frontal channel
m1<-lmer(lgFrontal~
           # Fixed-effects
           1+lgHz+
           # Random-effects
           (1+lgHz|subID:block), data=FILTERED, REML=FALSE)

Anova(m1, type="III")
summary(m1)
ranef(m1)
fixef(m1)
coef(m1)


dat<-data.frame(coef(m1)$'subID:block')
dat
write.csv(dat, "./data_EEG_initial_slopes_Frontal.csv")



## generating slopes for central channel
m1<-lmer(lgCentral~
           # Fixed-effects
           1+lgHz+
           # Random-effects
           (1+lgHz|subID:block), data=FILTERED, REML=FALSE)

Anova(m1, type="III")
summary(m1)
ranef(m1)
fixef(m1)
coef(m1)



dat<-data.frame(coef(m1)$'subID:block')
dat
write.csv(dat, "./data_EEG_initial_slopes_Central.csv")





## generating slopes for parietal channel
m1<-lmer(lgParietal~
           # Fixed-effects
           1+lgHz+
           # Random-effects
           (1+lgHz|subID:block), data=FILTERED, REML=FALSE)

Anova(m1, type="III")
summary(m1)
ranef(m1)
fixef(m1)
coef(m1)

dat<-data.frame(coef(m1)$'subID:block')
dat
write.csv(dat, "./data_EEG_initial_slopes_Parietal.csv")




## generating slopes for occipital channel
m1<-lmer(lgOccipital~
           # Fixed-effects
           1+lgHz+
           # Random-effects
           (1+lgHz|subID:block), data=FILTERED, REML=FALSE)

Anova(m1, type="III")
summary(m1)
ranef(m1)
fixef(m1)
coef(m1)

dat<-data.frame(coef(m1)$'subID:block')
dat
write.csv(dat, "./data_EEG_initial_slopes_Occipital.csv")

