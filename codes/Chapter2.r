#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


##############################################################
#Load the data
DF1 <- read.table("BaileyDensity.txt", header=TRUE)

#Subselect the data with MeanDepth > 800 and remove NAs
DF2 <- DF1[DF1$MeanDepth > 800,]
DF  <- na.exclude(DF2)

#Inspect the results
names(DF)
str(DF)

##############################################################
#House keeping: load packages and support functions
library(lattice)
library(nlme)
library(mgcv)
library(gstat)  #Trouble on a Mac (in 2015)

#Source our library file Highstatlib support file
source(file = "HighstatLibV8.R")
#This file can be downloaded from the book website.



##############################################################
#Section 2.2.First encounter with smoother

#Figure 2.1
DF$MyPch0               <- DF$Year
DF$MyPch0[DF$Year<1990] <- 1
DF$MyPch0[DF$Year>1990] <- 16

xyplot(Ykm ~ Xkm, 
       aspect = "iso", 
       col = 1, 
       data = DF,
       xlab = "X-coordinate", 
       ylab = "Y-coordinate",
       pch = DF$MyPch)


MyXLab <- "Mean sampling depth (m)"
MyYLab <- expression(paste("Fish density (",m^{-2}, ")"))

#Figure 2.2
par(mfrow=c(1,1), mar = c(5,5,2,3))
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab,
     #col= DF$Period,
     pch = DF$MyPch0,
     cex.lab = 1.5)

#Fit a linear regression model
M1 <- lm(Dens ~ MeanDepth, data = DF)
print(summary(M1), digits = 2, signif.stars = FALSE)


#Figure 2.3
E1 <- rstandard(M1)
F1 <- fitted(M1)
par(mfrow=c(1,2), mar = c(5,5,2,2))
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values", 
     ylab = "Residuals", 
     cex.lab = 1.5)
abline(h = 0, v = 0)

plot(x = DF$MeanDepth, 
     y = E1, 
     xlab = "Mean depth", 
     ylab = "Residuals", 
     cex.lab = 1.5)
abline(h = 0)


#Figure 2.4
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     cex.lab = 1.5,
     xlab = MyXLab ,
     ylab = MyYLab)
abline(M1,lwd = 5)


#Figure 2.5
M3 <- lm(Dens ~ poly(MeanDepth,3), data = DF)
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab)

D1 <- data.frame(MeanDepth=seq(from = min(DF$MeanDepth),
                               to = max(DF$MeanDepth),
                               length = 100))

P1 <- predict(M3, newdata = D1)
lines(x = D1$MeanDepth,   
      y = P1, 
      col = 1, 
      lwd = 5)


#Figure 2.6
M4 <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 5), 
          data = DF)

par(mar = c(5,5,3,3))
plot(M4, cex.lab = 1.5)


#Figure 2.7
MyPch <-vector(length = length(DF$MeanDepth))
MyPch[DF$MeanDepth < 2000] <- 1
MyPch[DF$MeanDepth > 3000] <- 1
MyPch[DF$MeanDepth >= 2000 & DF$MeanDepth <= 3000] <- 16

par(mfrow = c(1,2), mar = c(5.5,5,2,2))
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab,
     pch = MyPch,
     cex.lab = 1.5)

points(x = 2500,
       y = -0.001,
       pch = 17, 
       cex = 2)
x1 <- seq(0, 0.03, length = 50)
y1 <- rep(2000, length(x1))
lines(x = y1,
      y = x1,
      lwd = 3,
      lty = 2)
lines(x = y1 + 1000,
      y = x1,
      lwd = 3,
      lty = 2)

tapply(DF$Dens, MyPch,mean)
#Moving average
MyDepth <- seq(from = 800, to = 4850, length = 150)
MyMean  <- vector(length = 150)
k <- 1
for (i in MyDepth){
   MyPch <-vector(length = length(DF$MeanDepth))
   l1 <- i -500
   l2 <- i + 500
   MyPch[DF$MeanDepth < l1] <- 1
   MyPch[DF$MeanDepth > l2] <- 1
   MyPch[DF$MeanDepth >= l1 & DF$MeanDepth<=l2] <- 2
   m<-tapply(DF$Dens, MyPch, mean)
   MyMean[k] <- m[2]
   k <- k + 1
}

par(mar = c(5,5,2,2))
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab,
     pch = MyPch,
     cex.lab = 1.5)
lines(x = MyDepth,
      y = MyMean,
      lwd = 5)
#End of Figure 2.7



#Figure 2.8
L1 <- loess(Dens ~ MeanDepth, 
            data = DF,
            span = 0.75)

P1 <- predict(L1, se = TRUE)

plot(x = DF$MeanDepth, 
     y = DF$Dens,
     pch = 16,
     xlab = MyXLab,
     ylab = MyYLab)
lines(x = DF$MeanDepth, y = P1$fit, lwd = 3)
lines(x = DF$MeanDepth, y = P1$fit + 2 * P1$se.fit, lwd = 3, lty = 2)
lines(x = DF$MeanDepth, y = P1$fit - 2 * P1$se.fit, lwd = 3, lty = 2)
abline(h = 0)


#Figure 2.9
M1A <- loess(Dens ~ MeanDepth, data = DF, span = 0.1)
M1B <- loess(Dens ~ MeanDepth, data = DF, span = 0.5)
M1C <- loess(Dens ~ MeanDepth, data = DF, span = 1)

P1A <- predict(M1A, se = TRUE)
P1B <- predict(M1B, se = TRUE)
P1C <- predict(M1C, se = TRUE)

P   <- c(P1A$fit, P1B$fit, P1C$fit)
SEs <- c(P1A$se.fit, P1B$se.fit, P1C$se.fit)
MD  <- rep(DF$MeanDepth, 3)
ID  <- factor(rep(c("span = 0.1","span = 0.5","span = 1"), each = nrow(DF)))

Y <- rep(DF$Dens, 3)
X <- rep(DF$MeanDepth, 3)

#With polygons
xyplot(Y ~ X | ID, 
       col = 1,
       xlab = list(label = "Mean depth", cex = 1.5),
       ylab = list(label = "Density", cex = 1.5),
       layout = c(3,1),
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel = function(x, y, subscripts,...){
              y1  <- P[subscripts]
              x1  <- MD[subscripts]
              se1 <- SEs[subscripts]
              panel.polygon(c(x1, rev(x1)),
                 c(y1-2*se1, rev(y1+2*se1)),
                 col =grey(0.5),border=NULL,
                 density =50   )
              panel.lines(x1,y1, lwd = 3, col = 1)
              panel.points(x, y, col = 1)
              })


#End of Figure 2.9



#################################################################
#Section 2.3 Applying GAM with mgcv

#Fit the model
M4 <- gam(Dens ~ s(MeanDepth, fx =TRUE,k = 5), 
          data = DF)

summary(M4)
anova(M4)
plot(M4)


#Figure 2.10
MyData4 <- data.frame(MeanDepth = seq(from = 800,
                                      to = 4865,
                                      length = 100))
P4 <- predict(M4, 
              newdata = MyData4, 
              se = TRUE)

par(mar = c(5,5,3,3))
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab)
lines(x = MyData4$MeanDepth, y = P4$fit, lwd=3)
lines(x = MyData4$MeanDepth, y = P4$fit + 2 * P4$se.fit, lwd =3, lty = 2)
lines(x = MyData4$MeanDepth, y = P4$fit - 2 * P4$se.fit, lwd =3, lty = 2)


####################################################################
#Section 2.4 Cross validation

#Figure 2.11
M4A <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 3), data=DF)
M4B <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 4), data=DF)
M4C <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 5), data=DF)
M4D <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 6), data=DF)
M4E <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 8), data=DF)
M4F <- gam(Dens ~ s(MeanDepth, fx = TRUE, k = 10), data=DF)

P4A <- predict(M4A, newdata=MyData4, se=TRUE)
P4B <- predict(M4B, newdata=MyData4, se=TRUE)
P4C <- predict(M4C, newdata=MyData4, se=TRUE)
P4D <- predict(M4D, newdata=MyData4, se=TRUE)
P4E <- predict(M4E, newdata=MyData4, se=TRUE)
P4F <- predict(M4F, newdata=MyData4, se=TRUE)

AllP <- c(P4A$fit, P4B$fit,
          P4C$fit, P4D$fit,
          P4E$fit, P4F$fit)

AllX <- rep(MyData4$MeanDepth, 6)
AllSE <- c(P4A$se.fit, P4B$se.fit,
           P4C$se.fit, P4D$se.fit,
           P4E$se.fit, P4F$se.fit)

Names <- c("df = 2", "df = 3",
           "df = 4", "df = 5",
           "df = 7", "df = 9")
ID <- rep(Names, each = 100)

Y1 <- rep(DF$Dens, 6)
X1 <- rep(DF$MeanDepth, 6)
ID1 <- rep(Names, each = nrow(DF))

xyplot(Y1 ~ X1 | factor(ID1),
       xlab = list(label = "Mean depth", cex = 1.5),
       ylab = list(label = "Density", cex = 1.5),
       strip = strip.custom(bg = 'white',
                            par.strip.text = list(cex = 1.2)),
       scales = list(alternating = T,
                     x = list(relation = "same"),
                     y = list(relation = "same")),
       panel = function(x,y,subscripts,...){
             panel.points(x,y, col =1, pch = 1)
             ThsP <- ID1[subscripts][1]
             xx1 <- AllX[ID==ThsP]
             yy1 <- AllP[ID==ThsP]
             panel.lines(xx1, yy1, col = 1, lwd = 5)}
             )
#End of figure 2.11



#To apply cross validation
M5 <- gam(Dens ~ s(MeanDepth), data = DF)
summary(M5)


#Figure 2.12
par(mar = c(5,5,3,3))
plot(M5, cex.lab = 1.5)



#####################################################
# Section 2.5 Model validation
#Normality and Homogeneity
#Figure 2.13

E5 <- resid(M5)
F5 <- fitted(M5)

par(mfrow = c(1, 2), mar = c(5, 5, 3, 3))
hist(E5, 
     breaks = 25, 
     main = "", 
     xlab = "Residuals")
plot(x = F5, 
     y = E5, 
     xlab = "Fitted values", 
     ylab = "Residuals")
abline(h = 0, v = 0)


#Independence
#Figure 2.14
par(mfrow = c(1, 2), mar = c(5, 5, 3, 3))
plot(x = DF$MeanDepth, 
     y = E5, 
     xlab = "Depth", 
     ylab = "Residuals")
abline(h = 0)

boxplot(E5 ~ factor(DF$Period), 
        xlab = "Period", 
        ylab = "Residuals")
abline(h = 0)



#####################################################
#Figure 2.15
MyCex <- abs(E5) / max(abs(E5))
MyCol <- vector(length = length(E5))
MyCol[E5 > 0]  <- gray(0.5)
MyCol[E5 <= 0] <- gray(0.2)

xyplot(Ykm ~ Xkm,
       data = DF,
       main = list(label = "Residuals", cex = 1.5),
       xlab = list(label = "X-coordinates", cex = 1.5),
       ylab = list(label = "Y-coordinates", cex = 1.5),
       aspect = "iso",
       pch = 16,
       col = MyCol,
       cex = 3 * (MyCex)^(1/6))



#Figure 2.16
library(sp)
mydata <- data.frame(E5, DF$Xkm, DF$Ykm)
coordinates(mydata) <- c("DF.Xkm","DF.Ykm")
V5 <- variogram(E5 ~ 1, mydata)
plot(V5)



#Figure 2.17 Influential observation
M5 <- gam(Dens ~ s(MeanDepth), data = DF)
n  <- length(DF$Dens)
ID <- 1:n

MD <- data.frame(MeanDepth = seq(from = min(DF$MeanDepth),
                                 to = max(DF$MeanDepth), 
                                 length = 100))
P5 <- predict(M5, 
              newdata = MD, 
              type = "terms")
Dif <- vector(length = n)
for (i in 1:n) {
 print(i)
 M5.i <- gam(Dens ~ s(MeanDepth), data = DF, subset = ID[-i])
 P5.i <- predict(M5.i, newdata = MD, type = "terms")
 Dif[i]<-sum((P5[1:100]-P5.i[1:100])^2)
 }
 
# plot(Dif,
     # type = "h", 
     # xlab = "Observation number", 
     # ylab = "Change in smoother")

par(mar = c(5,6,3,3))
plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab,
     cex = 3*(Dif/max(Dif))^(1/4),
     pch = 16,
     cex.lab = 1.5)

# points(DF$MeanDepth[3], DF$Dens[3], col=2, pch=16)
# points(DF$MeanDepth[1], DF$Dens[1], col=2, pch=16)
# points(DF$MeanDepth[2], DF$Dens[2], col=2, pch=16)
# points(DF$MeanDepth[6], DF$Dens[6], col=2, pch=16)
#End of code Figure 2.17



###########################################################################
#Section 2.6 GAM with more covariates

DF$fPeriod <- factor(DF$Period)
M6 <- gam(Dens ~ s(MeanDepth) + fPeriod, data=DF)
summary(M6)


#Figure 2.18
par(mar = c(2,2,2,2))
vis.gam(M6, theta = 120, color = "heat")


#Figure 2.19
MD1 <- data.frame(MeanDepth = seq(from = 804, to = 4865, length = 100),
                  fPeriod = "1")

MD2 <- data.frame(MeanDepth = seq(from = 804, to = 4865, length = 100),
                  fPeriod = "2")

P1 <- predict(M6, newdata = MD1)
P2 <- predict(M6, newdata = MD2)

MyPch <- vector(length = length(DF$Dens))
MyPch[DF$Period == 1] <- 1
MyPch[DF$Period == 2] <-16

plot(x = DF$MeanDepth, 
     y = DF$Dens,
     xlab = MyXLab,
     ylab = MyYLab,
     pch = MyPch)
lines(x = MD1$MeanDepth, y = P1, lwd = 3)
lines(x = MD2$MeanDepth, y = P2, lwd = 3)


#Interactions
M7 <- gam(Dens ~ s(MeanDepth, by = fPeriod) + fPeriod,
          data = DF)
summary(M7)


#Figure 2.20
#For panels A and B
par(mfrow = c(2, 2))
plot(M7, select = 1, main = "Period 1")
plot(M7, select = 2, main = "Period 2")

#And plot the predicted values
MD7.1 <- data.frame(MeanDepth = seq(from = min(DF$MeanDepth),
                                    to = max(DF$MeanDepth), 
                                    length = 100),
                    fPeriod = "1")
MD7.2 <- data.frame(MeanDepth = seq(from = min(DF$MeanDepth),
                                    to = max(DF$MeanDepth), length=100),
                    fPeriod = "2")

P7.1 <- predict(M7, newdata=MD7.1)
P7.2 <- predict(M7, newdata=MD7.2)

#Panel C
plot(x = DF$MeanDepth[DF$Period == 1], 
     y = DF$Dens[DF$Period == 1],
     xlab = MyXLab,
     ylab = MyYLab,
     ylim = c(0, 0.03), 
     main = "Period 1")
lines(x = MD7.1$MeanDepth, y = P7.1, lwd=3)

#Panel D
plot(x = DF$MeanDepth[DF$Period == 2], 
     y = DF$Dens[DF$Period == 2],
     xlab = MyXLab,
     ylab = MyYLab,
     ylim = c(0, 0.03), 
     main = "Period 2")
lines(x = MD7.2$MeanDepth, y = P7.2, lwd = 3)
#End of code Figure 2.20


AIC(M6,M7)

#Interactions; second implementation
M8 <- gam(Dens ~ -1 +
                 s(MeanDepth, by = as.numeric(fPeriod == "1"))+
                 s(MeanDepth, by = as.numeric(fPeriod == "2")),
                  data = DF)
                  
                  
#Interactions; third implementation
M9 <- gam(Dens ~ s(MeanDepth)+
                 s(MeanDepth, by = as.numeric(fPeriod == "2")),
          data = DF)
summary(M9)


#Figure 2.21
par(mfrow = c(2,2), mar = c(5,6,3,2))
plot(M9, select = 1, main = "Period 1", cex.lab = 1.5)
plot(M9, select = 2, main = "Period 2", cex.lab = 1.5)

P9.1 <- predict(M9, newdata = MD7.1)
P9.2 <- predict(M9, newdata = MD7.2)

plot(x = DF$MeanDepth[DF$Period == 1], 
     y = DF$Dens[DF$Period == 1],
     xlab = MyXLab,
     ylab = MyXLab,
     ylim = c(0, 3.092395e-02), 
     main= "Period 1")
lines(x = MD7.1$MeanDepth, 
      y = P9.1, 
      lwd = 3)

plot(x = DF$MeanDepth[DF$Period == 2], 
     y = DF$Dens[DF$Period == 2],
     xlab = MyXLab,
     ylab = MyXLab,
     ylim = c(0, 3.092395e-02), 
     main = "Period 2")
lines(x = MD7.2$MeanDepth, y = P9.2, lwd = 3)
#End of Figure 2.21



#############################################################################
# Section 2.7 Transforming the density data

#Apply a square root transformation
DF$SQ.Dens <- sqrt(DF$Dens)
M10 <- gam(SQ.Dens ~ s(MeanDepth, by = fPeriod) + fPeriod,
                     data=DF)


#Figure 2.22
F10 <- fitted(M10)
E10 <- resid(M10)
plot(x = DF$MeanDepth,
     y = E10,
     xlab = "Mean depth", 
     ylab = "Residuals",
     pch = DF$MyPch)
abline(h = 0)


#######################################
#Section 2.8 Allowing for heterogeneity

IMD <- DF$MeanDepth
IMD[DF$MeanDepth < 2000] <- 1
IMD[DF$MeanDepth >= 2000] <- 2
IMD <- factor(IMD)

M11 <- gamm(Dens ~ s(MeanDepth, by = fPeriod) + fPeriod,
            weights = varIdent(form =~ 1 | IMD),
            data = DF)
summary(M11$lme)
summary(M11$gam)


#Figure 2.23
plot(M11$lme, col = 1, pch = 16, cex.lab = 1.5)


########################################################
#Section 2.9 Transforming and allowing for heterogeneity

#Figure 2.24
#log transformation
DF$L.Dens <- log(DF$Dens)
M1 <- gam(Dens ~ s(MeanDepth, by = factor(Period)) + factor(Period),
          data = DF)

M2 <- gam(SQ.Dens ~ s(MeanDepth, by = factor(Period)) + factor(Period),
          data = DF)

M3 <- gamm(SQ.Dens ~ s(MeanDepth, by = factor(Period)) + factor(Period),
           weights = varIdent(form =~ 1 | IMD),
           method = "REML",
           data = DF)

M4 <- gam(L.Dens ~ s(MeanDepth, by = factor(Period)) + factor(Period),
          data = DF)

#Predict from each model
#Depth values per period
MeanDepth.1 <- DF$MeanDepth[DF$Period == 1]
MeanDepth.2 <- DF$MeanDepth[DF$Period == 2]

ND1 <- data.frame(MeanDepth=seq(min(MeanDepth.1),max(MeanDepth.1),by=1),Period=1)
ND2 <- data.frame(MeanDepth=seq(min(MeanDepth.2),max(MeanDepth.2),by=1),Period=2)

PM1.1 <- predict(M1,newdata = ND1, se = TRUE)
PM1.2 <- predict(M1,newdata = ND2, se = TRUE)

PM2.1 <- predict(M2,newdata = ND1, se = TRUE)
PM2.2 <- predict(M2,newdata = ND2, se= TRUE)

PM3.1 <- predict(M3$gam,newdata = ND1, se = TRUE)
PM3.2 <- predict(M3$gam,newdata = ND2, se = TRUE)

PM4.1 <- predict(M4,newdata = ND1, se = TRUE)
PM4.2 <- predict(M4,newdata = ND2, se = TRUE)

AllY <- c(DF$Dens,
          DF$SQ.Dens,
          DF$SQ.Dens,
          DF$L.Dens)
AllX <- rep(DF$MeanDepth, 4)
MyNames<-c("Gaussian GAM",
           "Gaussian GAM, square root density",
           "Gaussian GAM, square root density, varIdent",
           "Gaussian GAM, log density")

AllID <- rep(MyNames, each = 147)
AllID <- factor(AllID, levels = MyNames)

AllPeriod <- rep(DF$Period,4)
AllPch <- AllPeriod
AllPch[AllPeriod == 1] <- 16
AllPch[AllPeriod == 2] <- 16

#x axis for the lines
AllMD.1 <- c(ND1$MeanDepth, ND1$MeanDepth,
             ND1$MeanDepth, ND1$MeanDepth)
AllMD.2 <- c(ND2$MeanDepth, ND2$MeanDepth,
             ND2$MeanDepth, ND2$MeanDepth)

#y axis for fitted lines
AllFit.1<- c(PM1.1$fit, PM2.1$fit,
             PM3.1$fit, PM4.1$fit)

AllFit.2<- c(PM1.2$fit, PM2.2$fit,
             PM3.2$fit, PM4.2$fit)
NewID.1 <- rep(MyNames,each = 4062)
NewID.2 <- rep(MyNames,each = 4033)

AllSE.1 <- c(PM1.1$se.fit, PM2.1$se.fit,
             PM3.1$se.fit, PM4.1$se.fit)

AllSE.2 <- c(PM1.2$se.fit, PM2.2$se.fit,
             PM3.2$se.fit, PM4.2$se.fit)

xyplot(AllY ~ AllX | AllID,
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "free")),
  xlab = list(label = "Depth (m)", cex = 1.5),
  ylab = list(label = "Fish density",cex = 1.5),
  panel=function(x,y, subscripts,...){
    ID <- AllID[subscripts][1]
    print(ID)
    x1 <- AllMD.1[NewID.1==ID]
    x2 <- AllMD.2[NewID.2==ID]
    y1 <- AllFit.1[NewID.1==ID]
    y2 <- AllFit.2[NewID.2==ID]
    se1<- AllSE.1[NewID.1==ID]
    se2<- AllSE.2[NewID.2==ID]
    panel.grid(h = -1, v = 2)

    #panel.lines(x1,y1+1.96*se1,lwd=3,lty=2,col=1)
    #panel.lines(x1,y1-1.96*se1,lwd=3,lty=2,col=1)

    #panel.lines(x2,y2+1.96*se2,lwd=3,lty=2,col=2)
    #panel.lines(x2,y2-1.96*se2,lwd=3,lty=2,col=2)
    z1.low <- y1 - 1.96*se1
    z1.up  <- y1 + 1.96*se1
    z2.low <- y2 - 1.96*se2
    z2.up  <- y2 + 1.96*se2

    panel.polygon(c(x1,rev(x1)),c(z1.low,rev(z1.up)),
                  col="grey65",border=NULL)
    panel.polygon(c(x2,rev(x2)),c(z2.low,rev(z2.up)),
                  col="indianred1", border=NULL)
    panel.points(x,y,col=1,   #AllPeriod[subscripts],
                 pch=AllPch[subscripts],cex=0.7)
    panel.lines(x1,y1,lwd=3,lty=1,col=1)
    panel.lines(x2,y2,lwd=3,lty=1,col=1)
    })
#End of Figure 2.24
