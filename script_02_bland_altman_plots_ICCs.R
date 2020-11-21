
library("ggplot2"); library("blandr"); library("ggpubr"); library("patchwork");
library ("ICC"); library("psych");

#### IIC correlations Frontal intercepts/ offset --------------------------------------
DATA1 <- read.csv("./data_EEG_COMB25_11182020.csv", 
                  header = TRUE,na.strings=c("","NA","na"))

colnames(DATA1)

# Make the FOOOF exponents negative to match the spectral slopes:
DATA1 <- DATA1 %>% mutate(slope_frontal = -slope_frontal,
                          slope_central = -slope_central,
                          slope_parietal = -slope_parietal,
                          slope_occipital = -slope_occipital)

head(DATA1)



REST<-subset(DATA1, time=="rest")
head(REST)

PRE <- subset(REST, block=="pretest")
POST <- subset(REST, block=="posttest")
head(PRE)

#### LMER Frontal channel spectral slopes

# Passes data to the blandr.statistics function to generate Bland-Altman statistics
statistics.results <- blandr.statistics(PRE$lgHz_central, POST$lgHz_central)

# Generates a ggplot, with title changed
g1 <- blandr.plot.ggplot( statistics.results , plotTitle = "Pre and Post-test Spectral slope: Central Region" ,
                          ciDisplay = FALSE , ciShading = FALSE ) + 
  ylim (-1.5,1.5) +
  labs(y = "Differences", x = "Means") 
g2 <-  g1 + theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        plot.caption=element_text(size=16, face = "bold"),  # caption
        axis.title.x=element_text(size=16, face = "bold"),  # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16),  # X axis text
        axis.text.y=element_text(size=16)) # Y axis text
print(g2)


### FOOOF slopes  frontal channel --------------------------------
# Passes data to the blandr.statistics function to generate Bland-Altman statistics
statistics.results <- blandr.statistics(PRE$slope_central, POST$slope_central)


# Generates a ggplot, with title changed
g3 <- blandr.plot.ggplot( statistics.results , plotTitle = "Pre and Post-test Aperiodic slope: Central Region" ,
                          ciDisplay = FALSE , ciShading = FALSE ) + 
                          ylim (-1.5,1.5) +
                          labs(y = "Differences", x = "Means") 
g4 <-  g3 + theme_bw() + 
  theme(plot.title=element_text(size=16, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        plot.caption=element_text(size=16),  # caption
        axis.title.x=element_text(size=16, face="bold"), # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16),  # X axis text
        axis.text.y=element_text(size=16)) # Y axis text
print(g4)


# combine the plots in patchwork:
g2/g4

 
# ICCs and Reliability Analyses ------------------------------------------------
## IIC correlations Frontal intercepts/ offset ---------------------------------

# Initial data visualziations: LMER Slopes
g1<-ggplot(REST, aes(x = subID, y = lgHz_frontal)) +
  geom_point(aes(fill=block), pch=21, size=2, stroke=1.25)
g2<-g1+scale_x_discrete(name = "Subject ID)") +
  scale_y_continuous(name = "LMER Slopes")

plot(g2)


# Initial data Visualizations: FOOOF Exponents 
g1<-ggplot(REST, aes(x = subID, y = slope_frontal)) +
  geom_point(aes(fill=block), pch=21, size=2, stroke=1.25)
g2<-g1+scale_x_discrete(name = "Subject ID)") +
  scale_y_continuous(name = "FOOOF Exponent")

plot(g2)




## ICCs for Resting Data -------------------------------------------------------
head(PRE)
head(POST)

# Making a loop to run the ICC for each outcome:
k=1

for (i in c(7:22)){
  print(paste("START", k))
  k<-k+1
  print(paste(colnames(PRE)[i], "~", colnames(POST)[i]))
  x<-data.frame(PRE[,i], POST[,i])
  print(ICC(x,missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE))
  print("STOP")
}




