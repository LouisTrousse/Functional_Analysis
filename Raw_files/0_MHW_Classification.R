####################################################################
#### Code accompanying:
#
# Gomez-Gras et al. 2021. Climate change transforms the functional identity in Mediterranean coralligenous assemblages. Ecology Letters

# Code to produce Figure 1, used to classify sites in MHW-impcted versus non-impacted sites.

# Script written by: Nathaniel Bensoussan and Daniel Gomez Gras

####################################################################

# set working directory

setwd("~/xxx")

#### Load required packages

library(ggplot2)
library(grid)
library(gridExtra)


### We load a  dataframe containing both the fitted 5th degree polynomial thermotolerance response curve (based on previous experimental data for both Paramuricea
#clavata and Corallium rubrum) and the observed MHWs over the monitored periods (in terms of duration and intensity). 

# * Details on how the thermotolerance curve used here was fitted can be found in the Materials and Methods section of the manuscript and in Appendix 1 / Fig. S3 from the Supplementary material. 

## ** Details on how duration (in days) and intensity (in degrees celsius) of MHWs were calculated for each assemblage can also be found on the Materials and Methods section of the manuscript. 

data <- read.csv2("MHW_classification.csv", sep=";", dec=",")

#### We produce a plot for every assemblage:

# Passe_cor 2003-2018

Passe_cor <- data[data$Site=="Passe_cor",]

p1 <- ggplot(Passe_cor, aes(x=Temp, y=MHWDays))+ 
  geom_line(aes(x=Passe_cor$Temp, y=Passe_cor$fitted), colour="black", linetype="longdash", size=0.8)+
  geom_point(aes(x=Passe_cor$Temp, y=Passe_cor$MHWDays, colour=factor(Impacted)), size=Passe_cor$Size*25)+
  geom_text(aes(x=Passe_cor$Temp-0.5, y=Passe_cor$MHWDays-0.5, label=Passe_cor$Year), colour="black")+
 scale_color_manual(values = c("orange","white"))+
   labs(x="T (C)", y="MHW Days")+ 
  ggtitle("Passe_cor") +
  theme(axis.line=element_line(size=.5,colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),
        panel.background=element_blank())+
  theme(axis.text.y = element_text(colour="black", size=10))+
  theme(axis.text.x = element_text(colour="black", size=10))+
  theme(axis.title.y=element_text(colour="black", size=11))+
  scale_y_continuous(limits =c(0,60), breaks= c(0,10,20,30,40,50,60),expand = c(0, 0))+
  scale_x_continuous(limits= c(22.5,28), breaks= c(23,24,25,26,27,28)) +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.title.x=element_text(colour="black", size=11))+
  theme(title=element_text(size=11.5),
        legend.text=element_text(size=12),
        strip.background=element_rect(colour="black",fill="white",size=2),
        strip.text=element_text(size=11.5,vjust=1))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())+
  coord_flip()+
  theme(legend.position="none")

##Pzzu_cor 2003 - 2018

Pzzu_cor <- data[data$Site=="Pzzu_cor",]

p2 <- ggplot(Pzzu_cor, aes(x=Temp, y=MHWDays))+ 
  geom_line(aes(x=Pzzu_cor$Temp, y=Pzzu_cor$fitted), colour="black", linetype="longdash", size=0.8)+
  geom_point(aes(x=Pzzu_cor$Temp, y=Pzzu_cor$MHWDays, colour=factor(Impacted)), size=Pzzu_cor$Size*25)+
  geom_text(aes(x=Pzzu_cor$Temp-0.5, y=Pzzu_cor$MHWDays-0.5, label=Pzzu_cor$Year), colour="black")+
  scale_color_manual(values = c("orange","red","white"))+
  labs(x="T (C)", y="MHW Days")+ 
  ggtitle("Pzzu_cor") +
  theme(axis.line=element_line(size=.5,colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),
        panel.background=element_blank())+
  theme(axis.text.x = element_text(colour="black", size=10))+
  theme(axis.title.y=element_blank())+
  theme(axis.text.y=element_blank())+
  scale_y_continuous(limits =c(0,60), breaks= c(0,10,20,30,40,50,60),expand = c(0, 0))+
  scale_x_continuous(limits= c(22.5,28), breaks= c(23,24,25,26,27,28)) +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.title.x=element_text(colour="black", size=11))+
  theme(title=element_text(size=11.5),
        legend.text=element_text(size=12),
        strip.background=element_rect(colour="black",fill="white",size=2),
        strip.text=element_text(size=11.5,vjust=1))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())+
  coord_flip()+
  theme(legend.position="none")


##Pzzinu_par 2006 - 2016

Pzzinu_par <- data[data$Site=="Pzzinu_par",]



p3 <- ggplot(Pzzinu_par, aes(x=Temp, y=MHWDays))+ 
  geom_line(aes(x=Pzzinu_par$Temp, y=Pzzinu_par$fitted), colour="black", linetype="longdash", size=0.8)+
  geom_point(aes(x=Pzzinu_par$Temp, y=Pzzinu_par$MHWDays, colour=factor(Impacted)), size=Pzzinu_par$Size*25)+
  geom_text(aes(x=Pzzinu_par$Temp-0.5, y=Pzzinu_par$MHWDays-0.5, label=Pzzinu_par$Year), colour="black")+
  scale_color_manual(values = c("orange","white"))+
  labs(x="T (C)", y="MHW Days")+ 
  ggtitle("Pzzinu_par") +
  theme(axis.line=element_line(size=.5,colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),
        panel.background=element_blank())+
  theme(axis.text.x = element_text(colour="black", size=10))+
  theme(axis.title.y=element_blank())+
  theme(axis.text.y=element_blank())+
  scale_y_continuous(limits =c(0,60), breaks= c(0,10,20,30,40,50,60),expand = c(0, 0))+
  scale_x_continuous(limits= c(22.5,28), breaks= c(23,24,25,26,27,28)) +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.title.x=element_text(colour="black", size=11))+
  theme(title=element_text(size=11.5),
        legend.text=element_text(size=12),
        strip.background=element_rect(colour="black",fill="white",size=2),
        strip.text=element_text(size=11.5,vjust=1))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())+
  coord_flip()+
  theme(legend.position="none")



##Pzzu_par 2006 - 2018

Pzzu_par <- data[data$Site=="Pzzu_par",]

p4 <- ggplot(Pzzu_par, aes(x=Temp, y=MHWDays))+ 
  geom_line(aes(x=Pzzu_par$Temp, y=Pzzu_par$fitted), colour="black", linetype="longdash", size=0.8)+
  geom_point(aes(x=Pzzu_par$Temp, y=Pzzu_par$MHWDays, colour=factor(Impacted)), size=Pzzu_par$Size*25)+
  geom_text(aes(x=Pzzu_par$Temp-0.5, y=Pzzu_par$MHWDays-0.5, label=Pzzu_par$Year), colour="black")+
  scale_color_manual(values = c("orange","red","white"))+
  labs(x="T (C)", y="MHW Days")+
  ggtitle("Pzzu_par") +
  theme(axis.line=element_line(size=.5,colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),
        panel.background=element_blank())+
  theme(axis.text.x = element_text(colour="black", size=10))+
  theme(axis.title.y=element_blank())+
  theme(axis.text.y=element_blank())+
  scale_y_continuous(limits =c(0,60), breaks= c(0,10,20,30,40,50,60),expand = c(0, 0))+
  scale_x_continuous(limits= c(22.5,28), breaks= c(23,24,25,26,27,28)) +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.title.x=element_text(colour="black", size=11))+
  theme(title=element_text(size=11.5),
        legend.text=element_text(size=12),
        strip.background=element_rect(colour="black",fill="white",size=2),
        strip.text=element_text(size=11.5,vjust=1))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())+
  coord_flip()+
  theme(legend.position="none")


###Gabin_par 1999 - 2009 


Gabin_par <- data[data$Site=="Gabin_par",]

p5 <- ggplot(Gabin_par, aes(x=Temp, y=MHWDays))+ 
  geom_line(aes(x=Gabin_par$Temp, y=Gabin_par$fitted), colour="black", linetype="longdash", size=0.8)+
  geom_point(aes(x=Gabin_par$Temp, y=Gabin_par$MHWDays, colour=factor(Impacted)), size=Gabin_par$Size*25)+
  geom_text(aes(x=Gabin_par$Temp+0.8, y=Gabin_par$MHWDays-0.5, label=Gabin_par$Year), colour="black")+
  scale_color_manual(values = c("orange","red","white"))+
  labs(x="T (C)", y="MHW Days")+ 
  ggtitle("Gabin_par") +
  #coord_cartesian(ylim =c(23,28))+
  theme(axis.line=element_line(size=.5,colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),
        panel.background=element_blank())+
  theme(axis.text.x = element_text(colour="black", size=10))+
  theme(axis.title.y=element_blank())+
  theme(axis.text.y=element_blank())+
  scale_y_continuous(limits =c(0,60), breaks= c(0,10,20,30,40,50,60),expand = c(0, 0))+
  scale_x_continuous(limits= c(22.5,28), breaks= c(23,24,25,26,27,28)) +
  #ggtitle("CC Direct effects") +
  theme(plot.title = element_text(hjust=0.5))+
  theme(axis.title.x=element_text(colour="black", size=11))+
  #geom_smooth(aes(colour = factor(Location)),method = "lm", se = FALSE, size=1.5)+
  theme(title=element_text(size=11.5),
        legend.text=element_text(size=12),
        strip.background=element_rect(colour="black",fill="white",size=2),
        strip.text=element_text(size=11.5,vjust=1))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())+
  coord_flip()+
  theme(legend.position="none")


##We make the combined plot to produce Figure 1.

pdf("Figure_1.pdf", width=13, height=2.8)
grid.arrange(p1,p2,p3,p4,p5, ncol=5,nrow=1, widths =  c(2.4,2,2,2,2))
dev.off() 

################# End of the code.
