library(lme4); library (ICC); library(psych); 
library(car); library(tidyverse); library(lmerTest)

# Multilevel models for LMER slopes versus FOOOF -------------------------------
DATA3 <- read.csv("./data_EEG_COMB1v1.csv", 
                  header = TRUE,na.strings=c("","NA","na"),
                  stringsAsFactors = TRUE)
head(DATA3)

# Make the FOOOF exponents negative to match the spectral slopes:
DATA3 <- DATA3 %>% mutate(slope_frontal = -slope_frontal,
                          slope_central = -slope_central,
                          slope_parietal = -slope_parietal,
                          slope_occipital = -slope_occipital)



PRAC_wide <-subset(DATA3, time == "prac")
head(PRAC_wide)

SLOPES_WIDE<-PRAC_wide %>% select(-c("intercept_frontal", "intercept_central", "intercept_parietal", 
                                  "intercept_occipital", "offset_frontal", "offset_central", 
                                  "offset_parietal", "offset_occipital"))
head(SLOPES_WIDE)

SLOPES_LONG <-gather(data = SLOPES_WIDE, key = method, value = slopes, c("lgHz_frontal", "lgHz_central", "lgHz_parietal", "lgHz_occipital",
                                                                         "slope_frontal", "slope_central", "slope_parietal", "slope_occipital"))

head(SLOPES_LONG)


INT_WIDE<-select(PRAC_wide, c("subID", "block", "intercept_frontal", "intercept_central", "intercept_parietal", 
                              "intercept_occipital", "offset_frontal", "offset_central", 
                              "offset_parietal", "offset_occipital"))
head(INT_WIDE)

INT_LONG <-gather(data = INT_WIDE, key = method, value = ints, c("intercept_frontal", "intercept_central", "intercept_parietal", 
                                                                   "intercept_occipital", "offset_frontal", "offset_central", 
                                                                   "offset_parietal", "offset_occipital"))

head(INT_LONG)

SLOPES_LONG$ints<-INT_LONG$ints
head(SLOPES_LONG)

# Splitting the "method" string into region and method variables:
SLOPES_LONG <- SLOPES_LONG %>% separate(col=method, 
                                       into = c("method","region"),
                                       sep = "_")
head(SLOPES_LONG)
summary(SLOPES_LONG$method)
SLOPES_LONG$method <- fct_recode(SLOPES_LONG$method, FOOOF = "slope", LMER = "lgHz")
SLOPES_LONG$method <- fct_relevel(SLOPES_LONG$method, "LMER", "FOOOF")

SLOPES_LONG <- SLOPES_LONG %>% arrange(subID, block)

write.csv(SLOPES_LONG, "./data_PRAC_long.csv")


data_PRAC_long1 <-subset(SLOPES_LONG, subID != "NA")


write.csv(data_PRAC_long1, "data_PRAC_long1.csv")


#### LMER model for PRAC -------------------------------
DATA <- read.csv("./data_PRAC_long1.csv", 
                 header = TRUE,na.strings=c("","NA","na"),
                 stringsAsFactors = TRUE)

head(DATA)


DATA$method.c<-(as.numeric(DATA$method)-1.5)*(-1) ## convert method to numerical

DATA$blocknum <-as.numeric(as.character(substr(DATA$block, 2, 3)))
summary(DATA$blocknum)

DATA$block.c<-DATA$blocknum-mean(DATA$blocknum)

head(DATA)

DATA$diff.c <- DATA$diff - mean(DATA$diff)


# Overall LMER Model -----------------------------------------------------------
M1 <- lmer(slopes ~ 1+ block.c*diff.c*method*region+
             (1|subID)+(1|subID:method)+(1|subID:region), 
           data = DATA, REML = FALSE)

Anova(M1, type="III")
summary(M1)

## Post Hoc Tests --------------------------------------------------------------
## Post hoc tests for the frontal channels -------------------------------------

FRONTAL <- subset(DATA, region == "frontal")
head(FRONTAL)


Fron_mod <- lmer(slopes ~ 1+ block.c*method+diff.c*method+
                   (1+block.c+diff.c+method|subID), 
                 data = FRONTAL, REML = FALSE,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=5e5)))

Anova(Fron_mod, type="III")
summary(Fron_mod)


Fron_mod_LMER <- lmer(slopes ~ 1+ block.c+diff.c+
                   (1+block.c+diff.c|subID), 
                 data = FRONTAL[FRONTAL$method=="LMER",], REML = FALSE,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=5e5)))
Anova(Fron_mod_LMER, type="III")
summary(Fron_mod_LMER)

Fron_mod_FOOOF <- lmer(slopes ~ 1+ block.c+diff.c+
                        (1+block.c+diff.c|subID), 
                      data = FRONTAL[FRONTAL$method=="FOOOF",], REML = FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
Anova(Fron_mod_FOOOF, type="III")
summary(Fron_mod_FOOOF)


# Post hoc tests for LMER method for central channels --------------------------
CENTRAL <- subset(DATA, region == "central")
head(CENTRAL)

Cent_mod <- lmer(slopes ~ 1+ block.c*method+diff.c*method+
                   (1+block.c+diff.c+method|subID), 
                 data = CENTRAL, REML = FALSE,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=5e5)))

Anova(Cent_mod, type="III")
summary(Cent_mod)


Cent_mod_LMER <- lmer(slopes ~ 1+ block.c+diff.c+
                        (1+block.c+diff.c|subID), 
                      data = CENTRAL[CENTRAL$method=="LMER",], REML = FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
Anova(Cent_mod_LMER, type="III")
summary(Cent_mod_LMER)


Cent_mod_FOOOF <- lmer(slopes ~ 1+ block.c+diff.c+
                        (1+block.c+diff.c|subID), 
                      data = CENTRAL[CENTRAL$method=="FOOOF",], REML = FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
Anova(Cent_mod_FOOOF, type="III")
summary(Cent_mod_FOOOF)



# Post hoc tests for LMER method for parietal channels --------------------------
PARIETAL <- subset(DATA, region == "parietal")
head(PARIETAL)

Par_mod <- lmer(slopes ~ 1+ block.c*method+diff.c*method+
                   (1+block.c+diff.c+method|subID), 
                 data = PARIETAL, REML = FALSE,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=5e5)))

Anova(Par_mod, type="III")
summary(Par_mod)



Par_mod_LMER <- lmer(slopes ~ 1+ block.c+diff.c+
                        (1+block.c+diff.c|subID), 
                      data = PARIETAL[PARIETAL$method=="LMER",], REML = FALSE,
                      control=lmerControl(optimizer="bobyqa",
                                          optCtrl=list(maxfun=5e5)))
Anova(Par_mod_LMER, type="III")
summary(Par_mod_LMER)


Par_mod_FOOOF <- lmer(slopes ~ 1+ block.c+diff.c+
                       (1+block.c+diff.c|subID), 
                     data = PARIETAL[PARIETAL$method=="FOOOF",], REML = FALSE,
                     control=lmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=5e5)))
Anova(Par_mod_FOOOF, type="III")
summary(Par_mod_FOOOF)


# Post hoc tests for LMER method for occipital channels --------------------------
OCCIPITAL <- subset(DATA, region == "occipital")
head(OCCIPITAL)

Occ_mod <- lmer(slopes ~ 1+ block.c*method+diff.c*method+
                  (1+block.c+diff.c+method|subID), 
                data = OCCIPITAL, REML = FALSE,
                control=lmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=5e5)))

Anova(Occ_mod, type="III")
summary(Occ_mod)



Occ_mod_LMER <- lmer(slopes ~ 1+ block.c+diff.c+
                       (1+block.c+diff.c|subID), 
                     data = OCCIPITAL[OCCIPITAL$method=="LMER",], REML = FALSE,
                     control=lmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=5e5)))
Anova(Occ_mod_LMER, type="III")
summary(Occ_mod_LMER)


Occ_mod_FOOOF <- lmer(slopes ~ 1+ block.c+diff.c+
                       (1+block.c+diff.c|subID), 
                     data = OCCIPITAL[OCCIPITAL$method=="FOOOF",], REML = FALSE,
                     control=lmerControl(optimizer="bobyqa",
                                         optCtrl=list(maxfun=5e5)))
Anova(Occ_mod_FOOOF, type="III")
summary(Occ_mod_FOOOF)



# Figure 5 ---------------------------------------------------------------------
head(DATA)
DATA$region <- fct_recode(DATA$region, Frontal = "frontal", Central = "central",
                            Parietal = "parietal", Occipital = "occipital")
DATA$region <- fct_relevel(DATA$region, "Frontal", "Central", "Parietal", "Occipital")


# Figure 5A: LMER Slopes by Block and Region
DAT2 <- DATA %>% filter(method=="LMER") %>%
              group_by(region, block) %>%
              summarize(blocknum = blocknum[1],
                        mean_slope = mean(slopes))
head(DAT2)

ggplot(DAT2, aes(x = blocknum, y = mean_slope)) +
  geom_line(aes(col=region, lty=region), lwd=1)+
  scale_x_continuous(name = "Block of Practice", limits=c(0,20)) +
  scale_y_continuous(name = "Mean Spectral Slope", limits=c(-1.8, -1.0)) +
  theme_bw() + labs(col="Region", lty="Region")+
  theme(axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = "bottom") 


# Figure 5B: LMER Slopes by Difficulty and Region
DAT3 <- DATA %>% filter(method=="LMER") %>%
  group_by(region, diff) %>%
  summarize(mean_slope = mean(slopes))

head(DAT3)

ggplot(DAT3, aes(x = diff, y = mean_slope)) +
  geom_line(aes(col=region, lty=region), lwd=1)+
  scale_x_continuous(name = "Difficulty of Practice", breaks=c(0:9)) +
  scale_y_continuous(name = "Mean Spectral Slope", limits=c(-1.8, -1.0)) +
  theme_bw() + labs(col="Region", lty="Region")+
  theme(axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = "bottom") 




# Figure 5C: FOOOF Slopes by Block and Region
DAT4 <- DATA %>% filter(method=="FOOOF") %>%
  group_by(region, block) %>%
  summarize(blocknum = blocknum[1],
            mean_slope = mean(slopes))
head(DAT4)

ggplot(DAT4, aes(x = blocknum, y = mean_slope)) +
  geom_line(aes(col=region, lty=region), lwd=1)+
  scale_x_continuous(name = "Block of Practice", limits=c(0,20)) +
  scale_y_continuous(name = "Mean Aperiodic Slope", limits=c(-1.8, -1.0)) +
  theme_bw() + labs(col="Region", lty="Region")+
  theme(axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = "bottom") 



# Figure 5D: FOOOF Slopes by Difficulty and Region
DAT5 <- DATA %>% filter(method=="FOOOF") %>%
  group_by(region, diff) %>%
  summarize(mean_slope = mean(slopes))

head(DAT5)

ggplot(DAT5, aes(x = diff, y = mean_slope)) +
  geom_line(aes(col=region, lty=region), lwd=1)+
  scale_x_continuous(name = "Difficulty of Practice", breaks=c(0:9)) +
  scale_y_continuous(name = "Mean Aperiodic Slope", limits=c(-1.8, -1.0)) +
  theme_bw() + labs(col="Region", lty="Region")+
  theme(axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = "bottom") 




# Analysis of the Agreement between Spectral and Aperiodic slopes during the Task
head(DATA)

DAT6 <- DATA %>% group_by(subID, method, region) %>%
  summarize(mean_slope = mean(slopes))

head(DAT6)

table(DAT6$method, DAT6$subID)

DAT7<-spread(data = DAT6, key = method, value=mean_slope)

ggplot(DAT7, aes(x = FOOOF, y = LMER)) +
  geom_point(aes(col=region), pch=21, size=2) +  
  stat_smooth(aes(col=region), method="lm", lwd=0.5)+
  facet_wrap(~region) +
  scale_x_continuous(name = "Aperiodic Slope (FOOOF)") +
  scale_y_continuous(name = "Spectral Slope (LMER)") +
  theme_bw() + labs(col="Region", lty="Region")+
  theme(axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = "none") 


ICC(DAT7[DAT7$region=="Frontal", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
ICC(DAT7[DAT7$region=="Central", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
ICC(DAT7[DAT7$region=="Parietal", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
ICC(DAT7[DAT7$region=="Occipital", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)




# ICCs between FOOOF and LMER for Resting data ----
DAT8 <- read.csv("./data_EEG_COMB1v1.csv", 
                  header = TRUE, na.strings=c("","NA","na"),
                 stringsAsFactors = TRUE)
head(DAT8)

# Make the FOOOF exponents negative to match the spectral slopes:
DAT8 <- DAT8 %>% mutate(slope_frontal = -slope_frontal,
                          slope_central = -slope_central,
                          slope_parietal = -slope_parietal,
                          slope_occipital = -slope_occipital)

head(DAT8)


REST <- DAT8 %>% filter(time=="rest") %>%
  select(subID, lgHz_frontal, slope_frontal, lgHz_central, slope_central,
         lgHz_parietal, slope_parietal, lgHz_occipital, slope_occipital) %>%
  group_by(subID) %>%
  summarize(lgHz_frontal=mean(lgHz_frontal), 
            slope_frontal=mean(slope_frontal), 
            lgHz_central=mean(lgHz_central), 
            slope_central=mean(slope_central),
            lgHz_parietal=mean(lgHz_parietal), 
            slope_parietal=mean(slope_parietal), 
            lgHz_occipital=mean(lgHz_occipital), 
            slope_occipital=mean(slope_occipital)
            )
head(REST)


REST_LONG <-gather(data = REST, key = method, value = slopes, c("lgHz_frontal", "lgHz_central", "lgHz_parietal", "lgHz_occipital",
                                                                         "slope_frontal", "slope_central", "slope_parietal", "slope_occipital"))

head(REST_LONG)


# Splitting the "method" string into region and method variables:
REST_LONG <- REST_LONG %>% separate(col=method, 
                                        into = c("method","region"),
                                        sep = "_")
head(REST_LONG)

REST_LONG$region <- fct_recode(REST_LONG$region, Frontal = "frontal", Central = "central",
                          Parietal = "parietal", Occipital = "occipital")
REST_LONG$region <- fct_relevel(REST_LONG$region, "Frontal", "Central", "Parietal", "Occipital")

REST_LONG$method <- fct_recode(REST_LONG$method, FOOOF = "slope", LMER = "lgHz")
REST_LONG$method <- fct_relevel(REST_LONG$method, "LMER", "FOOOF")

head(REST_LONG)

REST_LONG<-spread(data = REST_LONG, key = method, value=slopes)
head(REST_LONG)

ggplot(REST_LONG, aes(x = FOOOF, y = LMER)) +
  geom_point(aes(col=region), pch=21, size=2) +  
  stat_smooth(aes(col=region), method="lm", lwd=0.5)+
  facet_wrap(~region) +
  scale_x_continuous(name = "Aperiodic Slope (FOOOF)") +
  scale_y_continuous(name = "Spectral Slope (LMER)") +
  theme_bw() + labs(col="Region", lty="Region")+
  theme(axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=16, colour="black", face="bold"),
        strip.text.x = element_text(size = 16),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face="bold"),
        legend.position = "none") 


ICC(REST_LONG[REST_LONG$region=="Frontal", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
ICC(REST_LONG[REST_LONG$region=="Central", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
ICC(REST_LONG[REST_LONG$region=="Parietal", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)
ICC(REST_LONG[REST_LONG$region=="Occipital", c("FOOOF", "LMER")], missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE)


