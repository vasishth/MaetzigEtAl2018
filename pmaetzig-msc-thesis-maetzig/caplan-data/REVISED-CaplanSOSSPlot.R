library(reshape)
library(plyr)
library(ggplot2)

source('WuKaiserVasishthCogSci2016/R/magnifytext.R')
source('WuKaiserVasishthCogSci2016/R/multiplot.R')
source('WuKaiserVasishthCogSci2016/R/plotresults.R')
source('WuKaiserVasishthCogSci2016/R/regionmeans.R')
source('WuKaiserVasishthCogSci2016/R/regionplot.R')

# rewrote the 'regionmeans' function using reshape2 and plyr
source('regionmeans_new.R')

dat <- read.table('rcdata.txt', header = TRUE, dec = ',')

# basic sanity checks
dat.EC <- subset(dat, grp=="EC")
dat.Aph <- subset(dat, grp=="LCVA")
round(tapply(dat.EC$rt, INDEX=list(dat.EC$seg,dat.EC$sent), mean, na.rm=TRUE))
round(tapply(dat.Aph$rt, INDEX=list(dat.Aph$seg,dat.Aph$sent), mean, na.rm=TRUE))
xtabs(~subj+seg, dat.Aph)
xtabs(~subj+item+sent, dat.Aph)

# set coding of 'seg' to region type, get correct order of factor levels
dat$roi <- ifelse(dat$seg == "a", "NP1",
                  ifelse(dat$seg == "b", "who",
                         ifelse(dat$seg == "c", ifelse(dat$senttype == "SO", "NP2", "V1"),
                                ifelse(dat$seg == "d", ifelse(dat$senttype == "SO", "V1", "NP2"),
                                       ifelse(dat$seg == "e", "V2",
                                              ifelse(dat$seg == "f", "NP3", NA))))))
dat$roi<-as.factor(as.character(dat$roi))
dat$roi <- factor(dat$roi, levels=c("NP1", "who", "NP2", "V1", "V2", "NP3"))
dat$rt <- as.numeric(dat$rt)

# sanity check (no NA should appear as factor level)
levels(dat$roi)

# calculating region means using regionmeans_new function
dat.regionmeans.EC <- regionmeans_new(subset(dat, grp=="EC"))
dat.regionmeans.Aph <- regionmeans_new(subset(dat, grp!="EC"))

regionplot(dat.regionmeans.EC, maintitle="Caplan SO/SS, elderly controls", legendpos=c(0, 0.375))
regionplot(dat.regionmeans.Aph, maintitle="Caplan SO/SS, aphasics", legendpos=c(0, 0.375))

# removing the NP3 roi because of ridiculous SE
dat.regionmeans.EC.NoNP3 <- regionmeans_new(subset(dat, grp=="EC" & roi!="NP3"))
dat.regionmeans.Aph.NoNP3 <- regionmeans_new(subset(dat, grp!="EC" & roi!="NP3"))

# Plotting without NP3 roi
regionplot(dat.regionmeans.EC.NoNP3, maintitle="Caplan SO/SS, elderly controls", legendpos=c(0, 0.375))
regionplot(dat.regionmeans.Aph.NoNP3, maintitle="Caplan SO/SS, aphasics", legendpos=c(0, 0.375))
