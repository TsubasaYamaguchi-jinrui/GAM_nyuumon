#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


##############################################################
#Load the data
DF1 <-read.table("BaileyDensity.txt", 
                 header = TRUE)
DF2 <- DF1[DF1$MeanDepth > 800,]
DF  <- na.exclude(DF2)

#Inspect the results
names(DF)
str(DF)


##############################################################
#Housekeeping: load packages and support functions
library(lattice) #For multi-panel graphs
library(MASS)    #For NB GLM
library(mgcv)    #For GAMs



###############################################
#Section 4.2 review of GLM. Count data approach

DF$fPeriod <- factor(DF$Period)


#Figure 4.1
par(mfrow = c(2,2), mar = c(5,5,2,2))
dotchart(DF$SweptArea, 
         cex.lab = 1.5,
         xlab = "Swept area",
         ylab = "Order of the data")

dotchart(DF$TotAbund, 
         cex.lab = 1.5,
         xlab = "Total abundance",
         ylab = "Order of the data")

plot(x = DF$MeanDepth, 
     y = DF$SweptArea,
     xlab = "Mean depth (m)",
     ylab = "Swept area",
     cex.lab = 1.5)

plot(x = DF$MeanDepth, 
     y = DF$TotAbund,
     xlab = "Mean depth (m)",
     ylab = "Total abundance",
     cex.lab = 1.5)

######################################################
#Section 4.4 Results of Piosson and NB GLM

#Fit the models C1 and C2
C1 <- glm(TotAbund ~ MeanDepth * fPeriod,
          family = poisson, 
          data = DF)

C2 <- glm(TotAbund ~ MeanDepth * fPeriod,
          family = poisson(link = "identity"),
          start = coef(C1),
          data = DF)
#Warnings in C2!

AIC(C1, C2)

#Calculate the overdispersion
E1 <- resid(C1, type = "pearson")
N  <- nrow(DF)
p  <- length(coef(C1))
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion

MyXLab <- "Mean sampling depth (m)"
MyYLab <- expression(paste("Fish density (",m^{-2}, ")"))
DF$MyPch <- c(1,16)[DF$Period]
#Figure 4.3
par(mar = c(5,5,2,2))
plot(x = DF$MeanDepth, 
     y = E1,
     cex.lab = 1.5,
     pch = DF$MyPch,
     xlab = MyXLab,
     ylab = "Pearson residuals")
abline(h = 0)

#Fit model C3
C3 <- glm.nb(TotAbund ~ MeanDepth * fPeriod,
             data = DF)


####################################################################
#Section 4.6 Poisson GLM with offset

#Fit Model C4
DF$LogSA <- log(DF$SweptArea)
C4 <- glm(TotAbund ~ MeanDepth * fPeriod + offset(LogSA),
         family = poisson, 
         data = DF)


####################################################################
#Section 4.8 Poisson and NB GAM with offset
#Fit Model C6
C6 <- gam(TotAbund ~ s(MeanDepth, by = fPeriod) + fPeriod + offset(LogSA),
          family = poisson, 
          data = DF)

#Calcualte overdispersion
E6 <- resid(C6, type = "pearson")
Overdispersion <- sum(E6^2) / (C6$df.res)
Overdispersion

summary(C6)


#To plot the smoothers we can use the following commands
par(mfrow = c(1,2))
plot(C6)

par(mfrow = c(2,2))
plot(C6, all.terms = TRUE)

par(mfrow = c(1,2))
plot(C6,select = c(1), main ="Period 1")
plot(C6,select = c(2), main ="Period 2")

par(mfrow = c(2,2))
plot(C6, all.terms = TRUE, resid = TRUE, scale = FALSE)

par(mfrow = c(1,2))
plot(C6, shade = TRUE, scale = FALSE)


#Fit Model C7
C7 <- gam(TotAbund ~ s(MeanDepth, by = fPeriod) + fPeriod +
                     offset(LogSA),
          family = negbin(theta = c(1,10)),
          data = DF)
summary(C7)
#See the warning for the more up to date mgcv version.
#Use:
C7 <- gam(TotAbund ~ s(MeanDepth, by = fPeriod) + fPeriod +
                     offset(LogSA),
          family = nb(),
          data = DF)
#theta <- getTheta(C7)


#Fit Model C8
C8 <- gam(TotAbund ~ s(MeanDepth) + fPeriod +
                     offset(LogSA),
          family = nb(theta = 2.048),  #Slightly difference values as in book
          data = DF)                   #Due to different mgcv version! 

AIC(C7, C8)


#Figure 4.4
par(mfrow=c(1,2), mar = c(5,6,3,2))
plot(C7, select = 1, main = "Period 1")
plot(C7, select = 2, main = "Period 2")


##########################################################################
#Figure 4.5
MeanDepth.1 <- DF$MeanDepth[DF$fPeriod=="1"]
MeanDepth.2 <- DF$MeanDepth[DF$fPeriod=="2"]

MLSA1 <- mean(DF$LogSA[DF$fPeriod=="1"])
MLSA2 <- mean(DF$LogSA[DF$fPeriod=="2"])

ND1 <- data.frame(MeanDepth = seq(min(MeanDepth.1),
                                  max(MeanDepth.1),
                                  by = 1),
                  fPeriod = "1", 
                  LogSA = MLSA1)
ND2 <- data.frame(MeanDepth = seq(min(MeanDepth.2),
                                  max(MeanDepth.2),
                                  by = 1),
                  fPeriod = "2", 
                  LogSA = MLSA2)

PM1.1 <- predict(C7,newdata = ND1, se = TRUE)
PM1.2 <- predict(C7,newdata = ND2, se = TRUE)

AllY <- DF$TotAbun
AllX <- DF$MeanDepth

AllPeriod <- DF$Period
AllPch <- c(1,16)[AllPeriod]

#x axis for the lines
AllMD.1 <- ND1$MeanDepth
AllMD.2 <- ND2$MeanDepth

#y axis for fitted lines
AllFit.1 <- c(PM1.1$fit)
AllFit.2 <- c(PM1.2$fit)

AllSE.1 <- c(PM1.1$se.fit)
AllSE.2 <- c(PM1.2$se.fit)

xyplot(AllY ~ AllX ,
       strip = function(bg='white', ...) strip.default(bg='white', ...),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "free")),
       xlab = "Mean depth (m)",
       ylab = "Total abundance",
       panel=function(x,y, subscripts,...){
              x1 <- AllMD.1
              x2 <- AllMD.2
              y1 <- AllFit.1
              y2 <- AllFit.2
              se1<- AllSE.1
              se2<- AllSE.2
              panel.grid(h=-1, v= 2)

              z1.low <- exp(y1 - 1.96*se1)
              z1.up  <- exp(y1 + 1.96*se1)
              z2.low <- exp(y2 - 1.96*se2)
              z2.up  <- exp(y2 + 1.96*se2)

              panel.polygon(c(x1,rev(x1)),c(z1.low,rev(z1.up)),
                            col=grey(0.4),border=NULL)
              panel.polygon(c(x2,rev(x2)),c(z2.low,rev(z2.up)),
                            col=grey(0.7), border=NULL)
              panel.points(x,y, col = 1,
                           pch=AllPch[subscripts],cex=1)
              panel.lines(x1,exp(y1),lwd=3,lty=1,col=1)
              panel.lines(x2,exp(y2),lwd=3,lty=1,col=1)
            })
#I guess ggplot2 code would be much shorter!





#Figure 4.6
#Model validation should be for model M7
#Not M8!
E7 <- resid(C7, type = "pearson")
F7 <- fitted(C7, type = "response")

par(mfrow = c(2,2), mar = c(5,5,2,3))
plot(x = F7, 
     y = E7, 
     xlab = "Fitted values", 
     ylab = "Pearson Residuals")

plot(x = DF$MeanDepth, 
     y = E7, 
     xlab = "Mean depth (m)", 
     ylab ="Pearson residuals")
abline(h = 0, lty = 2)

plot(E7 ~ fPeriod, 
     data = DF, 
     ylab = "Pearson residuals", 
     xlab = "Period")
     
hist(E7, 
     xlab = "Pearson residuals", 
     main = "")

############## End of code #########