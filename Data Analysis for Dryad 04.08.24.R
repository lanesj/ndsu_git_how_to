rm(list=ls())

#Load required packages
library(Rmisc) 
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(lubridate)
library(glmm)
library(predictmeans)
library(lsmeans)
library(gridExtra)
library(effects)
library(emmeans)
library(broom)
library(modelr)
library(performance)
library(ggpubr)

##############################################################################################################
##############################################################################################################

#Read in full data set
telo.dat<-read.csv("rTLDryad.masterfile.csv")

#Turn the date column into a real date
telo.dat$date = as.Date(telo.dat$date, "%m/%d/%Y")

#Use lubridate to add day of year to data frame 
telo.dat$day=yday(telo.dat$date)

#Calculate z-scores from the rTL's
telo.dat$rTL.z<-((telo.dat$rTL)-(mean(telo.dat$rTL)))/sd(telo.dat$rTL)

#Check correlation between z-scores and raw data
ggscatter(telo.dat, x = "rTL", y = "rTL.z", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "rTL", ylab = "Z-scores")

#Categorize year
telo.dat$year<-as.factor(telo.dat$year)

#Categorize assay plate number
telo.dat$plate<-as.factor(telo.dat$plate)

#Mean center nestling age (days). 
telo.dat$mean.centered.age<-(telo.dat$Age..day.)-(mean(telo.dat$Age..day., na.rm=TRUE))

#Mean center brood size. 
telo.dat$mean.centered.brood<-(telo.dat$total.brood.size)-(mean(telo.dat$total.brood.size, na.rm=TRUE))

#######################################################################################################################
#######################################################################################################################
#################################---------Site and telomere comparison----------######################################
#######################################################################################################################
#######################################################################################################################

#Full model for comparison of adult and nestlings relative telomere lengths 
site.lm<-lmer(rTL.z~site+year+Age..year.+(1|NID)+(1|plate), data = telo.dat)

#Check model
check_model(site.lm)

hist(residuals(site.lm))

#overall effect
anova(site.lm)

#Parwise comparision 
summary(site.lm)

telo.dat = telo.dat %>% filter(site!="NA")

#Make summary table
site.z.sum<-summarySE(telo.dat, measurevar="rTL.z", groupvars=c("site" ,"Habitat"), na.rm = T)

#Make figure
sz<-ggplot(telo.dat, aes(x=site, y=rTL.z))+ 
  geom_jitter(position=position_jitter(w=0.1, h=0.1), shape=1, size=2,  alpha = 0.5 )+
  geom_errorbar(data = site.z.sum, mapping = aes(x = site, y = rTL.z,
                                                ymin = rTL.z - ci, ymax = rTL.z+ ci),
                linewidth =.5, width=.20,color='black') + 
  geom_point(data = site.z.sum, mapping = aes(x = site, y= rTL.z),color='black',
             size=2.5,shape=16) + 
  facet_grid(~Habitat, scales ="free")+
  theme_bw()+
  theme(text =element_text(size =13))+
  theme(panel.grid.major =element_blank(), panel.grid.minor =element_blank())+
  theme(strip.placement = "outside")+
  theme(text =element_text(size =13))+
  ylab("Relative Telomere Length (z-scores)") +
  xlab("Site")

sz + theme(axis.title.x = element_text(size=15),
           axis.text.x  = element_text(size=12),
           axis.text.y  = element_text(size=12),
           axis.title.y = element_text(size=15))

#######################################################################################################################
#######################################################################################################################
#################################---------Adult and Nestling comparison----------######################################
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
#################################---------Adult and nestling comparison----------######################################
#######################################################################################################################

#Full model for comparison of adult and nestlings relative telomere lengths 
hy.hab<-lmer(rTL.z~Habitat*Age..year.+year+(1|NID)+(1|plate), data = telo.dat)

#Check model
check_model(hy.hab)

#summary output
summary(hy.hab)

#Interaction is not significant so dropping the interaction an rerun model. 
hy.hab<-lmer(rTL.z~Habitat+Age..year.+year+(1|NID)+(1|plate), data = telo.dat)

#Check model
check_model(hy.hab)

#summary output
summary(hy.hab)

#Habitat not significant so dropping variable. 
hy.hab<-lmer(rTL.z~Age..year.+year+(1|NID)+(1|plate), data = telo.dat)

#Check model
check_model(hy.hab)

#r2 values
r2(hy.hab)

#Conditional R2: 0.465
#Marginal R2: 0.075

#summary output
summary(hy.hab)

#################################---------Adult and nestling comparison figure:z scores----------######################################

#Make summary table
age.z.sum<-summarySE(telo.dat, measurevar="rTL.z", groupvars=c("Age..year."), na.rm = T)

#Change variable titles to be more informative in both Data frames. 
telo.dat$Age..year.<- factor(telo.dat$Age..year., levels=c("HY","ASY"))

levels(telo.dat$Age..year.) <- c("Nestling","Adult")

age.z.sum$Age..year. <- factor(age.z.sum$Age..year., levels=c("HY","ASY"))

levels(age.z.sum$Age..year.) <- c("Nestling","Adult")

#Make figure
bz<-ggplot(telo.dat, aes(x=Age..year., y=rTL.z))+ 
  geom_jitter(position=position_jitter(w=0.1, h=0.1), shape=1, size=2,  alpha = 0.5 )+
  geom_errorbar(data = age.z.sum, mapping = aes(x = Age..year., y = rTL.z,
                                              ymin = rTL.z - se, ymax = rTL.z+ se),
                linewidth =.5, width=.20,color='black') + 
  geom_point(data = age.z.sum, mapping = aes(x = Age..year., y= rTL.z),color='black',
             size=2.5,shape=16) + 
  theme_classic()+
  theme(text =element_text(size =13))+
  ylab("Relative Telomere Length (z-scores)") +
  xlab("Age")

bz + theme(axis.title.x = element_text(size=15),
           axis.text.x  = element_text(size=12),
           axis.text.y  = element_text(size=12),
           axis.title.y = element_text(size=15))

############################################################################################################################
############################################################################################################################
######################---------Adult Data Analysis----------################################################################
############################################################################################################################
############################################################################################################################

#Filter to just adults
asy.dat<-filter(telo.dat, Age..year.=="Adult")

############################################################################################################################
#################################----------Adult difference between habitats----------######################################
############################################################################################################################

#run full model examining variation in adult relative telomere length between habitat types and sexes. 
asy.hab<-lmer(rTL.z~Habitat*sex+year+(1|plate), data = asy.dat)

#Check mdoel
check_model(asy.hab)

#summary output
summary(asy.hab)

#Drop interaction as it is not signifcant 
asy.hab<-lmer(rTL.z~Habitat+sex+year+(1|plate), data = asy.dat)

#Check assumptions of model
check_model(asy.hab)

#r2 values
r2(asy.hab)

#Conditional R2: 0.176
#Marginal R2: 0.050

#summary output
summary(asy.hab)

#################################---------Adult difference in rTL's between habitats figure:z-scores----------##########################################

#Summarize variables.
asy.z.sum<-summarySE(asy.dat, measurevar="rTL.z", groupvars=c("Habitat","sex"), na.rm = TRUE)

#Change variable titles and organize to be more informative in both data frames. 
asy.dat$sex<- factor(asy.dat$sex, levels=c("female","male"))

levels(asy.dat$sex) <- c("Female","Male")

asy.z.sum$sex <- factor(asy.z.sum$sex, levels=c("female","male"))

levels(asy.z.sum$sex) <- c("Female","Male")

#create a dodge for the graph
pd <- position_dodge(1)

#Create figure
az=ggplot(data=asy.dat, aes(y=rTL.z, x=sex, color = Habitat), na.rm = TRUE)+ 
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, jitter.height=0.1, dodge.width = 1), shape=1, size=2.5, alpha=.5)+
  geom_errorbar(data = asy.z.sum, mapping = aes(x = sex, y = rTL.z,
                                              ymin =rTL.z - se, ymax = rTL.z + se),
                linewidth=.70, width=.3, position = pd) +
  geom_point(data = asy.z.sum, mapping = aes(x =sex, y= rTL.z),
             size=2, position = pd) + 
  xlab("Sex")+
  ylab("Relative Telomere Length (z-scores)")+
  theme_classic()+
  labs(color = "Habitat") +
  scale_color_manual(values = c("Urban" = "gray60", "Rural" = "gray1"))

az + theme(axis.title.x = element_text(size=15),
          axis.text.x  = element_text(size=12),
          axis.text.y  = element_text(size=12),
          axis.title.y = element_text(size=15),legend.position="none")

#################################---------Adult and nestlings, main text----------######################################

#Filter to nests that are not parasatized
hy.asy.all<-filter(telo.dat, nest.para.fig!="yes")

#Make summary table
age.hab.sum<-summarySE(hy.asy.all, measurevar="rTL", groupvars=c("Habitat","Age..year."), na.rm = T)

#Make figure
bb<-ggplot(hy.asy.all, aes(x=Habitat, y=rTL, group = Age..year.))+ 
  geom_jitter(position=position_jitter(w=0.1, h=0.1), shape=1, size=2,  alpha = 0.5 )+
  geom_errorbar(data = age.hab.sum, mapping = aes(x = Habitat, y = rTL,
                                                  ymin = rTL - se, ymax = rTL+ se),
                linewidth =.5, width=.20,color='black') + 
  geom_point(data = age.hab.sum, mapping = aes(x = Habitat, y= rTL),color='black',
             size=2.5,shape=16) + 
  facet_grid(~Age..year., scales ="free")+
  theme_bw()+
  theme(text =element_text(size =13))+
  theme(panel.grid.major =element_blank(), panel.grid.minor =element_blank())+
  theme(strip.placement = "outside")+
  ylab("Relative Telomere Length (z-scores)") 


bb + theme(axis.title.x = element_text(size=15),
           axis.text.x  = element_text(size=12),
           axis.text.y  = element_text(size=12),
           axis.title.y = element_text(size=15))

##############################################################################################################
##############################################################################################################
##########################---------Nestlings Data Analysis----------##########################################
##############################################################################################################
##############################################################################################################

#Filter to just nestlings
hy.dat<-filter(telo.dat, Age..year.=="Nestling")

##############################################################################################################
###########---------Differences in nestling rTL between Habitat Types----------################################
##############################################################################################################

#Filter to nests that are not parasatized
hy.hab.dat<-filter(hy.dat, nest.para=="no")

#run full lmm for habitat type effects on nestling TL
hy.hab<-lmer(rTL.z~Habitat+mean.centered.age+mean.centered.brood+year+(1|NID)+(1|plate), data = hy.hab.dat)

#Check assumptions of model
check_model(hy.hab)

#r2 values
r2(hy.hab)

#Conditional R2: 0.694
#Marginal R2: 0.286

#summary output
summary(hy.hab)

#################################---------Nestling rTL by Habitat Type Figure: z-scores----------##########################################

#Summarize variables.
hy.z.sum<-summarySE(hy.hab.dat, measurevar="rTL.z", groupvars=c("Habitat"), na.rm = TRUE)

#Change variable titles to be more informative in both Data frames. 
hy.hab.dat$Habitat <- factor(hy.hab.dat$Habitat, levels=c("Rural","Urban"))

levels(hy.hab.dat$Habitat) <- c("Rural; Not Parasitized","Urban; Not Parasitized")

hy.z.sum$Habitat <- factor(hy.z.sum$Habitat, levels=c("Rural","Urban"))

levels(hy.z.sum$Habitat) <- c("Rural; Not Parasitized","Urban; Not Parasitized")

#Make figure
cz<-ggplot(hy.hab.dat, aes(x=Habitat, y=rTL.z))+ 
  geom_jitter(position=position_jitter(w=0.1, h=0.1), shape=1, size=2,  alpha = 0.5 )+
  geom_errorbar(data = hy.z.sum, mapping = aes(x = Habitat, y = rTL.z,
                                             ymin = rTL.z - se, ymax = rTL.z+ se),
                linewidth=.5, width=.20,color='black') + 
  geom_point(data = hy.z.sum, mapping = aes(x = Habitat, y= rTL.z),color='black',
             size=2.5,shape=16) + 
  theme_classic()+
  theme(text =element_text(size =13))+
  ylab("Relative Telomere Length (z-scores)") +
  xlab("Habitat Type")

cz + theme(axis.title.x = element_text(size=15),
           axis.text.x  = element_text(size=12),
           axis.text.y  = element_text(size=12),
           axis.title.y = element_text(size=15))

##############################################################################################################
##############---------The effect of brood parasitism on nestling rTL's----------#############################
##############################################################################################################

#Filter out nests that do not contain brood parasitism
hy.para.dat<-filter(hy.dat, para.comp=="yes")

#run lmm for brod para on TL 
hy.para<-lmer(rTL.z~bhco.hab+mean.centered.age+mean.centered.brood+year+site+(1|NID)+(1|plate), data =hy.para.dat)

#Check model
check_model(hy.para)

#summary output
summary(hy.para)

#Drop nestling age as it is not significant
hy.para<-lmer(rTL.z~bhco.hab+mean.centered.age+year+site+(1|NID)+(1|plate), data =hy.para.dat)

#Check model
check_model(hy.para)

#summary output
summary(hy.para)

#Drop Brood size
hy.para<-lmer(rTL.z~bhco.hab+year+site+(1|NID)+(1|plate), data =hy.para.dat)

#Check model
check_model(hy.para)

#r2 values
r2(hy.para)

#Conditional R2: 0.616
#Marginal R2: 0.152

#summary output
summary(hy.para)

#################################---------Nestling rTL's Brood parasitism Figure: z-scores----------##########################################

#Summarize variables.
hy.z.sum<-summarySE(hy.para.dat, measurevar="rTL.z", groupvars=c("bhco.hab"), na.rm = TRUE)

#Change variable titles to be more informative in both Data frames. 
hy.para.dat$bhco.hab <- factor(hy.para.dat$bhco.hab, levels=c("urban.no","urban.yes"))

levels(hy.para.dat$bhco.hab) <- c("Urban; Not Parasitized","Urban; Parasitized")

hy.z.sum$bhco.hab <- factor(hy.z.sum$bhco.hab, levels=c("urban.no","urban.yes"))

levels(hy.z.sum$bhco.hab) <- c("Urban; Not Parasitized","Urban; Parasitized")

#Make figure
zd<-ggplot(hy.para.dat, aes(x=bhco.hab, y=rTL.z))+ 
  geom_jitter(position=position_jitter(w=0.1, h=0.1), shape=1, size=2,  alpha = 0.5 )+
  geom_errorbar(data = hy.z.sum, mapping = aes(x = bhco.hab, y = rTL.z,
                                             ymin = rTL.z - se, ymax = rTL.z+ se),
                linewidth=.5, width=.20,color='black') + 
  geom_point(data = hy.z.sum, mapping = aes(x = bhco.hab, y= rTL.z),color='black',
             size=2.5,shape=16) + 
  theme_classic()+
  theme(text =element_text(size =13))+
  ylab("Relative Telomere Length (z-scores)") +
  xlab("Parasitism Classification")

zd + theme(axis.title.x = element_text(size=15),
           axis.text.x  = element_text(size=12),
           axis.text.y  = element_text(size=12),
           axis.title.y = element_text(size=15))

