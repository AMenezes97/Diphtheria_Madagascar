## useful packages ## (you should click on package installer, etc to get them)
require(dplyr)
require(beeswarm)
require(RColorBrewer)



############################################################################
## Vacc calculation for the introduction, with 94% for PT and 97% for DT
expected.both <- 0.94*0.97; round(expected.both,3)
round(0.97*(1-0.94),3)
round(0.94*(1-0.97),3)
round(1-0.94*0.97-0.97*(1-0.94)-0.94*(1-0.97),3)
############################################################################

## 0. Bring in country level data and plot timeseries for Figure 1   ##################################################################################################

## time series data (from WHO)
dfts <- read.csv("diphtheria_june9.csv")
yr <- as.numeric(substring(colnames(dfts)[2:ncol(dfts)],2,5))

## vaccination data (likewise)
dfvacc <- read.csv("diphtheria.vacc-june9.csv")
yr.vacc <- as.numeric(substring(colnames(dfvacc)[2:ncol(dfvacc)],2,5))

## find index of countries of interest in case numbers
mada <- dfts[dfts[,1]=="Madagascar",2:ncol(dfts)]
india <- dfts[dfts[,1]=="India",2:ncol(dfts)]
indonesia <- dfts[dfts[,1]=="Indonesia",2:ncol(dfts)]

## pull out Madagascar vacc
mada.vacc <- dfvacc[dfvacc[,1]=="Madagascar",2:ncol(dfvacc)]

## plot cases in different settings for comparison
plot(yr,mada, type="b", xlab="", ylab="Reported diphtheria cases", ylim=range(c(mada, india, indonesia),na.rm=TRUE),pch=19)
points(yr,india, type="b", xlab="", ylab="Reported diphtheria cases",col="forestgreen",pch=19)
points(yr,indonesia, type="b", xlab="", ylab="Reported diphtheria cases", col="coral",pch=19)

## Plot out Figure 1
par(mfrow=c(1,1), bty="l",mar=c(5,5,4,5))
plot(yr,mada, type="b", xlab="", ylab="Reported diphtheria cases", pch=19, ylim=c(-2,3000),xlim=c(1980,2020),axes=FALSE)
polygon(c(2016-c(15,15,8),2016-c(8,15)),c(-2,3000,3000,-2,-2), col="grey", border="NA")
polygon(c(2016+c(-0.5,-0.5,0.5),2016+c(0.5,-0.5)),c(-2,3000,3000,-2,-2), col="grey", border="NA")
points(yr,mada, type="b",pch=19)
axis(1); axis(2)

points(yr.vacc,3000*(mada.vacc/100), type="b", xlab="", pch=19, col=4,xlim=c(1980,2020))
points(yr.vacc,3000*(mada.vacc/100), type="b", pch=19, col=4)
axis(4,at=c(0,1500,3000),lab=c(0,0.5,1),col=4,col.lab=4)
mtext("Vaccination coverage",4,line=2, col=4)


## 1. Bring in serology data and tidy it up  ##################################################################################################
df1 <- read.csv("Madagascar_Serology.csv")


## makes ages in years and months
df1$Age.complete <- df1$Age..months./12 + df1$Age..years.

## make categorical variable for DT where 0 means seronegative
m.anti.DT <- rep(NA,nrow(df1)) 
m.anti.DT[df1$anti.DT.levels == "<0.01 "] <- 0
m.anti.DT[df1$anti.DT.levels == "0.01-0.099"] <- 1
m.anti.DT[df1$anti.DT.levels == "0.1-1.0"] <- 2
m.anti.DT[df1$anti.DT.levels == ">1.0"] <- 3
df1$m.anti.DT <- as.ordered(m.anti.DT)

## make categorical variable for PT where <5 means seronegative
df1$anti.PT.cat<-cut(df1$anti.PT.Titer..IU.ml., br=c(-1,4.99,201)) 
m.anti.PT <- rep(NA,nrow(df1)) 
m.anti.PT[df1$anti.PT.cat == "(-1,4.99]"] <- 0
m.anti.PT[df1$anti.PT.cat == "(4.99,201]"] <- 1
df1$m.anti.PT <- as.ordered(m.anti.PT)


## tidy up some of the other variables
df1$Nutritionnal.status[df1$Nutritionnal.status=="" | df1$Nutritionnal.status=="Outlier"] <- NA
df1$District <- as.factor(df1$District)
df1$Nutritionnal.status <- as.factor(df1$Nutritionnal.status)


## function to put things on a 0-1 scale
inv.logit <- function(x) return(exp(x)/(1+exp(x)))

## make colors for plotting nutritional status
cols <- rep("black",nrow(df1))
cols[df1$Nutritionnal.status=="NA"] <- "grey"
cols[df1$Nutritionnal.status=="MM"] <- "blue"
cols[df1$Nutritionnal.status=="MS"] <- "red"

## make colors for districts
cols.dist <- rep("black",nrow(df1))
cols.dist[df1$District=="ANT"] <- brewer.pal(8,'Set2')[1]
cols.dist[df1$District=="ATS"] <- brewer.pal(8,'Set2')[2]
cols.dist[df1$District=="MAH"] <- brewer.pal(8,'Set2')[3]
cols.dist[df1$District=="MID"] <- brewer.pal(8,'Set2')[4]
cols.dist[df1$District=="TOL"] <- brewer.pal(8,'Set2')[5]

## make colors for DP only and both DP And PT
cols.vacc <- rep("black",nrow(df1))
cols.vacc[df1$m.anti.PT>0 & df1$m.anti.DT>0 ] <- "orange"
cols[df1$Nutritionnal.status=="MS"] <- "red"


## 2. focus on 6-11 month old ##################################################################################################
## 2.1 identify 0 dose and determinants 
#The following doesn't work with the updated code
fit1 <- glm(I(m.anti.DT==0)~ District, data = df1, subset=df1$Ages.group=="6-11 m")
summary(fit1)

##CIs on OR
exp(fit1$coeff[2]+c(-1,1)*1.96*summary(fit1)$coeff[2,2])
exp(fit1$coeff[3]+c(-1,1)*1.96*summary(fit1)$coeff[3,2])



beeswarm(log(pmax(df1$Anti.DT..IU.ml.[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL"],0.008))~
           as.character(df1$District[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL"]),
         axes=FALSE,xlab="District", ylab="Antibody titre (log scale)", pch=19, cex=0.7, xlim=c(0.5,4),
         pwcol=cols[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL"])
axis(side=1, at=c(1,2,3), lab=c("ANT","ATS","MID")); axis(side=2,  at=log(c(0.01,0.1,1)),lab=c(0.01,0.1,1))
abline(h=log(c(0.01,0.1,1)),lty=3, col=c("blue","grey","grey"))
text(c(1,2,3)+0.45,rep(log(0.008),3),round(inv.logit(fit1$coeff[1]+c(0,fit1$coeff[2:3])),2),cex=0.8) #needs to get removed




## 1.2 characterize heterogeneity in those vaccinated
##2.2 Characterize heterogeneity in those seropositive 

# group nutrit status (remove district, no effect)
fit2 <- lm(log(df1$Anti.DT..IU.ml.)~I(Nutritionnal.status=="NN"), data = df1, subset=df1$Ages.group=="6-11 m"& m.anti.DT>0)
summary(fit2)

round(fit2$coeff[1]+c(-1,1)*1.96*summary(fit2)$coeff[1,2],2)
round(fit2$coeff[2]+c(-1,1)*1.96*summary(fit2)$coeff[2,2],2)



beeswarm(log(pmax(df1$Anti.DT..IU.ml.[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL" & m.anti.DT>0]))~
           substring(as.character(df1$Nutritionnal.status[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL" & m.anti.DT>0]),1,1),
         axes=FALSE,xlab="Nutritional Status", ylab="Antibody titre (log scale)", pch=19, cex=0.7, xlim=c(0.5,3),
         pwcol=cols.dist[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL" & m.anti.DT>0])
#axis(side=1, at=c(1,2,3), lab=c("MM","MS","NN")); axis(side=2,  at=log(c(0.01,0.1,1)),lab=c(0.01,0.1,1))
axis(side=1, at=c(1,2), lab=c("Malnourished","Normal")); axis(side=2,  at=log(c(0.01,0.1,1)),lab=c(0.01,0.1,1))
abline(h=log(c(0.01,0.1,1)),lty=3, col=c("blue","grey","grey"))
points(c(1,2),c(fit2$coeff[1]+c(0,fit2$coeff[2])),pch=15, cex=2)

# group nutrit status (remove district, no effect) RESTRICT ANALYSIS TO INDIVIDUALS SEROPOSITIVE TO BOTH aka individuals likely to have been vaccinated 
fit2.1 <- lm(log(df1$Anti.DT..IU.ml.)~I(Nutritionnal.status=="NN"), data = df1, subset=df1$Ages.group=="6-11 m"& m.anti.DT>0 & m.anti.PT>0)
summary(fit2.1)

round(fit2.1$coeff[1]+c(-1,1)*1.96*summary(fit2.1)$coeff[1,2],2)
round(fit2.1$coeff[2]+c(-1,1)*1.96*summary(fit2.1)$coeff[2,2],2)

beeswarm(log(pmax(df1$Anti.DT..IU.ml.[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL" & m.anti.DT>0 & m.anti.PT>0]))~
           substring(as.character(df1$Nutritionnal.status[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL" & m.anti.DT>0 & m.anti.PT>0]),1,1),
         axes=FALSE,xlab="Nutritional Status", ylab="Antibody titre (log scale)", pch=19, cex=0.7, xlim=c(0.5,3),
         pwcol=cols.dist[df1$Ages.group=="6-11 m" & df1$District!="MAH" & df1$District!="TOL" & m.anti.DT>0 & m.anti.PT>0])
#axis(side=1, at=c(1,2,3), lab=c("MM","MS","NN")); axis(side=2,  at=log(c(0.01,0.1,1)),lab=c(0.01,0.1,1))
axis(side=1, at=c(1,2), lab=c("Malnourished","Normal")); axis(side=2,  at=log(c(0.01,0.1,1)),lab=c(0.01,0.1,1))
abline(h=log(c(0.01,0.1,1)),lty=3, col=c("blue","grey","grey"))
points(c(1,2),c(fit2.1$coeff[1]+c(0,fit2.1$coeff[2])),pch=15, cex=2)


## check on effects imm card, nutritional status...
fit3 <- lm(log(df1$Anti.DT..IU.ml.)~District+Nutritionnal.status+as.factor(Immunization.card), data = df1, subset=df1$Ages.group=="6-11 m"& m.anti.DT>0)
summary(fit3)


## check on effects imm card, nutritional status... RESTRICT ANALYSIS TO INDIVIDUALS SEROPOSITIVE TO BOTH aka individuals likely to have been vaccinated 
fit3.1 <- lm(log(df1$Anti.DT..IU.ml.)~District+Nutritionnal.status+as.factor(Immunization.card), data = df1, subset=df1$Ages.group=="6-11 m"& m.anti.DT>0 & m.anti.PT>0)
summary(fit3.1)



## 3. focus on 8-15 monthsyo ##################################################################################################

fit1a <- glm(I(m.anti.DT==0)~ District, data = df1, subset=df1$Ages.group=="8-15 y")
summary(fit1a)

exp(fit1a$coeff[1])/(1+exp(fit1a$coeff[1]))

exp(fit1a$coeff[2:5])
exp(fit1a$coeff[2:5]+1.96*summary(fit1a)$coeff[2:5,2])
exp(fit1a$coeff[2:5]-1.96*summary(fit1a)$coeff[2:5,2])


chs <- df1$Ages.group=="8-15 y"
beeswarm(log(pmax(df1$Anti.DT..IU.ml.[chs],0.008))~
           as.character(df1$District[chs]),
         axes=FALSE,xlab="District", ylab="Antibody titre (log scale)", pch=19, cex=0.7, xlim=c(0.5,5.5),
         pwcol=cols[chs])
axis(side=1, at=c(1:5), lab=c("ANT","ATS","MAH","MID","TOL")); axis(side=2,  at=log(c(0.01,0.1,1)),lab=c(0.01,0.1,1))
abline(h=log(c(0.01,0.1,1)),lty=3, col=c("blue","grey","grey"))


round(inv.logit(fit1a$coeff[1]+c(0,fit1a$coeff[2:5])),2)


fit2a <- lm(log(pmax(df1$Anti.DT..IU.ml.,0.008)) ~ District, data = df1, subset=df1$Ages.group=="8-15 y" & m.anti.DT>0)
summary(fit2a)


(fit2a$coeff[1])
(fit2a$coeff[1]+c(-1,1)*1.96*summary(fit2a)$coeff[1,2])

(fit2a$coeff[2:5])
(fit2a$coeff[2:5]+1.96*summary(fit2a)$coeff[2:5,2])
(fit2a$coeff[2:5]-1.96*summary(fit2a)$coeff[2:5,2])



## note that immunization card and age are no goes
fit3a <- lm(log(pmax(df1$Anti.DT..IU.ml.,0.008)) ~ Age..years. + District+as.factor(Immunization.card), data = df1, subset=df1$Ages.group=="8-15 y")
summary(fit3a)





## 4. Look at vaccination status among the different combinations. ############################################################################
#Look at seropositive status to diphtheria and pertussis 

## make categorical variable for PT where <5 means seronegative
df1$anti.PT.cat<-cut(df1$anti.PT.Titer..IU.ml., br=c(-1,4.99,201)) 
m.anti.PT <- rep(NA,nrow(df1)) 
m.anti.PT[df1$anti.PT.cat == "(-1,4.99]"] <- 0
m.anti.PT[df1$anti.PT.cat == "(4.99,201]"] <- 1

df1$m.anti.PT <- as.ordered(m.anti.PT)



## index to restrict to just youngest kids
chs <- df1$Ages.group=="6-11 m" & !is.na(df1$Ages.group)

## calculate and print out for the totals
tab1 <- table(df1$m.anti.PT[chs],1*(df1$m.anti.DT[chs]>0))
tab1
tab1/sum(tab1)

## calculate and print out numbers for Antananarivo
tab1.ANT <- table(df1$m.anti.PT[chs & df1$District=="ANT"],1*(df1$m.anti.DT[chs & df1$District=="ANT"]>0))
tab1.ANT

## calculate and print out numbers for ATS
tab1.ATS <- table(df1$m.anti.PT[chs & df1$District=="ATS"],1*(df1$m.anti.DT[chs & df1$District=="ATS"]>0))
tab1.ATS

## calculate and print out numbers for MID
tab1.MID <- table(df1$m.anti.PT[chs & df1$District=="MID"],1*(df1$m.anti.DT[chs & df1$District=="MID"]>0))
tab1.MID 




## index to restrict to just older individuals
ads <- df1$Ages.group=="8-15 y" & !is.na(df1$Ages.group)

## calculate and print out for the totals
tab2 <- table(df1$m.anti.PT[ads],1*(df1$m.anti.DT[ads]>0))
tab2
tab2/sum(tab2)

## calculate and print out numbers for Antananarivo
tab2.ANT <- table(df1$m.anti.PT[ads & df1$District=="ANT"],1*(df1$m.anti.DT[ads & df1$District=="ANT"]>0))
tab2.ANT

## calculate and print out numbers for ATS
tab2.ATS <- table(df1$m.anti.PT[ads & df1$District=="ATS"],1*(df1$m.anti.DT[ads & df1$District=="ATS"]>0))
tab2.ATS

## calculate and print out numbers for MID
tab2.MID <- table(df1$m.anti.PT[ads & df1$District=="MID"],1*(df1$m.anti.DT[ads & df1$District=="MID"]>0))
tab2.MID 

## calculate and print out numbers for MAH 
tab2.MAH <- table(df1$m.anti.PT[ads & df1$District=="MAH"],1*(df1$m.anti.DT[ads & df1$District=="MAH"]>0))
tab2.MAH 

## calculate and print out numbers for TOL 
tab2.TOL <- table(df1$m.anti.PT[ads & df1$District=="TOL"],1*(df1$m.anti.DT[ads & df1$District=="TOL"]>0))
tab2.TOL




## 5. Bring in district level vaccination coverage ##################################################################################################

## Reported administrative vaccination coverage data
df.vacc <- read.csv("Routine_Vaccination_2015-2017.csv")

## plot progress through time #light is past, dark recent, also shown by arrow, years are 2015, 2016, 2017
cols <- colorRampPalette(c(rgb(0,0,1,0), rgb(0,0,1,1)), alpha = TRUE)(6)[2:4] 

## plot out the five that we care about ######################
focal.districts <- c("Antananarivo Renivohitra","Antsalova","Mahajanga I","Midongy Atsimo","Toliara I")
chs <- match(focal.districts,df.vacc[,2])
#df.vacc[chs,2]
adj <- 0.1


lower.vacc.1 <- apply(df.vacc[chs,c(2,5,8)+1][,],1,min,na.rm=TRUE)
upper.vacc.1 <- apply(df.vacc[chs,c(2,5,8)+1][,],1,max,na.rm=TRUE)
med.vacc.1 <- apply(df.vacc[chs,c(2,5,8)+1][,],1,median,na.rm=TRUE)

lower.vacc.2 <- apply(df.vacc[chs,c(3,6,9)+1][,],1,min,na.rm=TRUE)
upper.vacc.2 <- apply(df.vacc[chs,c(3,6,9)+1][,],1,max,na.rm=TRUE)
med.vacc.2 <- apply(df.vacc[chs,c(3,6,9)+1][,],1,median,na.rm=TRUE)

lower.vacc.3 <- apply(df.vacc[chs,c(4,7,10)+1][,],1,min,na.rm=TRUE)
upper.vacc.3 <- apply(df.vacc[chs,c(4,7,10)+1][,],1,max,na.rm=TRUE)
med.vacc.3 <- apply(df.vacc[chs,c(4,7,10)+1][,],1,median,na.rm=TRUE)

par(mfrow=c(1,1),mar=c(10,5,2,2))

plot(1:5,med.vacc.1, ylab="Percent vaccinated", col=cols[1], xlab="", axes=FALSE,pch=19, ylim=c(0.4,1), xlim=c(0.5,5.5))
for (j in 1:5) points(c(j,j),c(lower.vacc.1[j],upper.vacc.1[j]), type="l", col=cols[1])

points(1:5+adj,med.vacc.2,pch=19, col=cols[2])
for (j in 1:5) points(c(j,j)+adj,c(lower.vacc.2[j],upper.vacc.2[j]), type="l", col=cols[2])

points(1:5+adj+adj,med.vacc.3,pch=19, col=cols[3])
for (j in 1:5) points(c(j,j)+adj+adj,c(lower.vacc.3[j],upper.vacc.3[j]), type="l", col=cols[3])
axis(2); axis(1,at=1:5, lab=df.vacc[chs,2], las=2, cex.axis=0.65)



## plot against fraction that are seropositive to both in each location against average vaccination coverage for third dose##

prop.both <- c(tab1.ANT[2,2]/sum(tab1.ANT),tab1.ATS[2,2]/sum(tab1.ATS),tab1.MID[2,2]/sum(tab1.MID))
n.tot <- c(sum(tab1.ANT),sum(tab1.ATS),sum(tab1.MID))
ci.low.prop.both <- prop.both-1.96*sqrt(prop.both*(1-prop.both)/n.tot)
ci.high.prop.both <- prop.both+1.96*sqrt(prop.both*(1-prop.both)/n.tot)

lower.vacc <- apply(df.vacc[chs,c(4,7,10)+1][c(1,2,4),],1,min)
upper.vacc <- apply(df.vacc[chs,c(4,7,10)+1][c(1,2,4),],1,max)
med.vacc <- apply(df.vacc[chs,c(4,7,10)+1][c(1,2,4),],1,median)

plot(med.vacc,prop.both,xlim=c(0.0,1.0),ylim=c(0,0.6),
     xlab="Percent vaccinated 2015-2017, third dose", 
     ylab="Proportion with antibodies to PT and DP", pch=19, )
for (j in 1:3) points(rep(med.vacc[j],2),c(ci.low.prop.both[j],ci.high.prop.both[j]),type="l",lty=3)
for (j in 1:3) points(c(lower.vacc[j],upper.vacc[j]),rep(prop.both[j],2),type="l",lty=3)


#points(med.vacc,med.vacc*expected.both, pch=19,col="red",ylim=c(0,0.8))
abline(0,1,col="red",lty=7)
points(c(0,0.5,1),expected.both*c(0,0.5,1),col="blue", type="l",lwd=2)


## plot against fraction that are DP only in each location against average vaccination coverage for third dose in 2015-2017##
prop.Donly <- c(tab1.ANT[1,2]/sum(tab1.ANT),tab1.ATS[1,2]/sum(tab1.ATS),tab1.MID[1,2]/sum(tab1.MID))
n.tot <- c(sum(tab1.ANT),sum(tab1.ATS),sum(tab1.MID))
ci.low.prop.Donly <- prop.Donly-1.96*sqrt(prop.Donly*(1-prop.Donly)/n.tot)
ci.high.prop.Donly <- prop.Donly+1.96*sqrt(prop.Donly*(1-prop.Donly)/n.tot)

lower.vacc <- apply(df.vacc[chs,c(4,7,10)+1][c(1,2,4),],1,min)
upper.vacc <- apply(df.vacc[chs,c(4,7,10)+1][c(1,2,4),],1,max)
med.vacc <- apply(df.vacc[chs,c(4,7,10)+1][c(1,2,4),],1,median)

plot(med.vacc,prop.Donly,xlim=c(0.0,1.0),ylim=c(0,0.6),
     xlab="Percent vaccinated 2015-2017, third dose", 
     ylab="Proportion with antibodies to DP only", pch=19, )
for (j in 1:3) points(rep(med.vacc[j],2),c(ci.low.prop.Donly[j],ci.high.prop.Donly[j]),type="l",lty=3)
for (j in 1:3) points(c(lower.vacc[j],upper.vacc[j]),rep(prop.Donly[j],2),type="l",lty=3)



##Raise pertussis seropositive cut off to 40 IU/mL
##Supplementary Table 3 and 4
## make categorical variable for PT where <40 means seronegative
df1$anti.PT.cat<-cut(df1$anti.PT.Titer..IU.ml., br=c(-1,39,201)) 
m.anti.PT <- rep(NA,nrow(df1)) 
m.anti.PT[df1$anti.PT.cat == "(-1,39]"] <- 0
m.anti.PT[df1$anti.PT.cat == "(39,201]"] <- 1

df1$m.anti.PT <- as.ordered(m.anti.PT)



chs <- df1$Ages.group=="6-11 m" & !is.na(df1$Ages.group)

## redo tables
## calculate and print out for the totals
tab1 <- table(df1$m.anti.PT[chs],1*(df1$m.anti.DT[chs]>0))
tab1
tab1/sum(tab1)

## calculate and print out numbers for Antananarivo
tab1.ANT <- table(df1$m.anti.PT[chs & df1$District=="ANT"],1*(df1$m.anti.DT[chs & df1$District=="ANT"]>0))
tab1.ANT

## calculate and print out numbers for ATS
tab1.ATS <- table(df1$m.anti.PT[chs & df1$District=="ATS"],1*(df1$m.anti.DT[chs & df1$District=="ATS"]>0))
tab1.ATS

## calculate and print out numbers for MID
tab1.MID <- table(df1$m.anti.PT[chs & df1$District=="MID"],1*(df1$m.anti.DT[chs & df1$District=="MID"]>0))
tab1.MID 



## index to restrict to just older individuals
ads <- df1$Ages.group=="8-15 y" & !is.na(df1$Ages.group)

## redo tables
## calculate and print out for the totals
tab2 <- table(df1$m.anti.PT[ads],1*(df1$m.anti.DT[ads]>0))
tab2
tab2/sum(tab2)

## calculate and print out numbers for Antananarivo
tab2.ANT <- table(df1$m.anti.PT[ads & df1$District=="ANT"],1*(df1$m.anti.DT[ads & df1$District=="ANT"]>0))
tab2.ANT

## calculate and print out numbers for ATS
tab2.ATS <- table(df1$m.anti.PT[ads & df1$District=="ATS"],1*(df1$m.anti.DT[ads & df1$District=="ATS"]>0))
tab2.ATS

## calculate and print out numbers for MID
tab2.MID <- table(df1$m.anti.PT[ads & df1$District=="MID"],1*(df1$m.anti.DT[ads & df1$District=="MID"]>0))
tab2.MID 

## calculate and print out numbers for MAH 
tab2.MAH <- table(df1$m.anti.PT[ads & df1$District=="MAH"],1*(df1$m.anti.DT[ads & df1$District=="MAH"]>0))
tab2.MAH 
