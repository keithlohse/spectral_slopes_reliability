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
# Switchin to Log Base 10 to match FOOOF Figures
head(DATA_AVE)
FIG1B <- DATA_AVE
FIG1B$lgHz <- log(FIG1B$Hz, 10)
FIG1B$ave_lg_central <- log(FIG1B$ave_central, 10)
FIG1B <- subset(FIG1B, block=="pretest")
FIG1B_no_alpha <- subset(FIG1B, Hz<8 | Hz>=13)
summary(lm(ave_lg_central~lgHz, data=FIG1B_no_alpha))

ggplot(data = FIG1B, 
       mapping = aes(x = lgHz, y = ave_lg_central)) +
  scale_x_continuous(name = "Log Frequency (log(Hz))") +
  scale_y_continuous(name = "Log Mean Central Power (log(uV^2))",
                     breaks=c(-0.6,-0.4,-0.2, 0.0, 0.2, 0.4, 0.6, 0.8))+
  annotate("rect", xmin=log(8,10), xmax=log(13, 10), ymin=-0.6, ymax=0.7, 
           alpha=0.2, fill="grey", col="black")+
  geom_line(col="black", lwd=1) +
  stat_smooth(data=FIG1B_no_alpha, method="lm", col="red", lty=2, se=FALSE)+
  theme_classic()+
  theme(axis.text=element_text(size=14, colour="black"), 
        strip.text=element_text(size=14, colour="black"), 
        axis.title=element_text(size=14,face="bold"),
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
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill="grey90",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))





#### Generating lm slopes and intercepts ----------------------------------
# Generating R^2 and coefficients for each participant to assess model Fit
# We will create an empty data frame to store our output in...
head(FILTERED)

FILTERED$Factor <- factor(paste(FILTERED$subID, FILTERED$time, FILTERED$block, sep="_"))
FILTERED$Factor <- factor(FILTERED$Factor)
head(FILTERED)
FILTERED <- subset(FILTERED, subID != "a25") # Removing missing participant

LIST <- as.vector(unique(FILTERED$Factor))

index<-c(1:length(LIST))
DAT2<-data.frame(index) 
head(DAT2)

for (i in 1:length(LIST)) {
  SAMP <- FILTERED[FILTERED$Factor==LIST[i],]
  
  DAT2$cluster[i]<-LIST[i]
  DAT2$subID[i] <- as.vector(SAMP$subID[i])
  DAT2$time[i] <- as.vector(SAMP$time[i])
  DAT2$block[i] <- as.vector(SAMP$block[i])
  # Frontal
  DAT2$frontal_intercept[i]<-lm(lgFrontal~lgHz, data=SAMP)$coefficients[1]
  DAT2$frontal_slope[i]<-lm(lgFrontal~lgHz, data=SAMP)$coefficients[2]
  DAT2$frontal_rsquared[i]<-summary(lm(lgFrontal~lgHz, data=SAMP))$adj.r.squared
  
  # Central
  DAT2$central_intercept[i]<-lm(lgCentral~lgHz, data=SAMP)$coefficients[1]
  DAT2$central_slope[i]<-lm(lgCentral~lgHz, data=SAMP)$coefficients[2]
  DAT2$central_rsquared[i]<-summary(lm(lgCentral~lgHz, data=SAMP))$adj.r.squared
  
  # Frontal
  DAT2$parietal_intercept[i]<-lm(lgParietal~lgHz, data=SAMP)$coefficients[1]
  DAT2$parietal_slope[i]<-lm(lgParietal~lgHz, data=SAMP)$coefficients[2]
  DAT2$parietal_rsquared[i]<-summary(lm(lgParietal~lgHz, data=SAMP))$adj.r.squared
  
  # Frontal
  DAT2$occipital_intercept[i]<-lm(lgOccipital~lgHz, data=SAMP)$coefficients[1]
  DAT2$occipital_slope[i]<-lm(lgOccipital~lgHz, data=SAMP)$coefficients[2]
  DAT2$occipital_rsquared[i]<-summary(lm(lgOccipital~lgHz, data=SAMP))$adj.r.squared
  
  print(LIST[i])
}

head(DAT2)
write.csv(DAT2, "./data_LMER_params.csv")
