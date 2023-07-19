#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.



##############################################################
#Set the working directory (on a Mac) and load the data
setwd("/Users/Highstat/applicat/HighlandStatistics/Courses/FilesOnlineFollowUpRegressionGLMGAM/FinalExercises")
HVS <- read.table(file = "HVS.txt", 
                  header = TRUE, 
                  dec = ".")

#Inspect the results
str(HVS)
names(HVS)
##############################################################



##############################################################
#Load packages and library files
library(lattice)  #Needed for multi-panel graphs
library(mgcv)     #Needed for smoothing curves in scatterplots

#source the file: HighstatLibV8.R
source(file = "/Users/Highstat/applicat/HighlandStatistics/Courses/Data/HighstatLibV8.R")
#This file can be dowloaded from the book website.
##############################################################



##############################################################
#Housekeeping
HVS$OrbitalV    <- HVS$MeanOrbitalVolume
HVS$fPopulation <- factor(HVS$Population)
HVS$LatAbs      <- HVS$AbsoluteLatitude
HVS$CC          <- HVS$CranialCapacity
HVS$FM          <- HVS$FMarea_intercondyle
HVS$Illuminance <- HVS$Minimum_Illuminance
HVS$Temperature <- HVS$Minimum_Temperature_celsius
HVS$fGender     <- factor(HVS$Gender)
##############################################################


##############################################################
#Section 1.4 Data exploration

colSums(is.na(HVS))
#There is one missing value. Drop it
HVS2 <- na.exclude(HVS)


#Outliers
#Figure 1.2
MyVar <- c("OrbitalV", "LatAbs", "CC", "FM",
           "Illuminance", "Temperature")
Mydotplot(HVS2[,MyVar])
##############################################################


##############################################################
#Collinearity
#Figure 1.3
MyVar <- c("LatAbs", "CC", "FM", "Illuminance", "Temperature")
Mypairs(HVS2[, MyVar])
##############################################################



##############################################################
#Multipanel boxplots
#Figure 1.4
MyVar2 <- c("LatAbs", "CC", "FM", "Illuminance", "Temperature")
Mybwplot(HVS2, MyVar2, "fGender")
##############################################################



##############################################################
#Multipanel scatterplots
#Figure 1.5
Myxyplot(HVS2, MyVar2, "OrbitalV")
##############################################################



##############################################################
#Section 1.5 Multiple linear regression

#Bivariate linear regression model
M1 <- lm(OrbitalV ~ LatAbs, data = HVS)
print(summary(M1), digits = 2, signif.stars = FALSE)

#Model matrix
model.matrix(M1)
model.matrix(~ LatAbs, data = HVS)

#Figure 1.6
par(mfrow = c(2,2), mar = c(5,5,3,3))
plot(x = HVS$LatAbs, 
     y = HVS$OrbitalV,
     xlab = "LatAbs",
     ylab = "OrbitalV",
     cex.lab = 1.5)

plot(x = HVS$LatAbs, 
     y = HVS$OrbitalV,
     xlab = "LatAbs",
     ylab = "OrbitalV",
     cex.lab = 1.5)
abline(M1)


plot(x = HVS$LatAbs, 
     y = HVS$OrbitalV,
     xlab = "LatAbs",
     ylab = "OrbitalV",
     cex.lab = 1.5)
abline(M1)

N <- nrow(HVS)
F1 <- fitted(M1)
for (i in 1:N){
	segments(HVS$LatAbs[i], HVS$OrbitalV[i],
	        HVS$LatAbs[i], F1[i],
	        lty = 2)
}
##############################################################



##############################################################
#Multiplying and calculating the inverse
X <- model.matrix(M1)
betaHat <- solve(t(X) %*% X) %*% t(X) %*%
           HVS$OrbitalV
betaHat
   
#Calculate residuals
e    <- resid(M1)
estd <- rstandard(M1)

#Calcualte the H matrix
H <- X %*% solve(t(X) %*% X) %*% t(X) 
hat(X)  #Not sure where I got hatmatrix function from
##############################################################



##############################################################
#leverage values for each onservation
#Figure 1.7
plot(diag(H), 
     type = "h",
     xlab = "Crania",
     ylab = "Leverage",
     cex.lab = 1.5)
##############################################################


##############################################################
#Cook distance
#Figure 1.8
plot(cooks.distance(M1), type = "h")

#or

plot (M1, which = 4)
##############################################################



##############################################################
#Standard errors
SE <- summary(M1)$sigma *
        sqrt(diag(solve(t(X) %*% X)))
SE

#Obtain confidence intervals
Z <- rbind(coef(M1) - 1.96 * SE,
           coef(M1) + 1.96 * SE)
rownames(Z) <- c("Lower bound","Upper bound")
Z

qt(1-0.05/2,df =55-2)

#Obtain the covariance matrix
covBeta <- vcov(M1)

#Create the data frame
MyData <- data.frame(LatAbs = seq(from = 0.02, 
                                  to = 65,
                                  length = 10))
X <- model.matrix(~LatAbs, data = MyData)

#Calculate the SEs
MyData$P <- X %*% coef(M1)
MyData$SEy.ci <- sqrt(diag(X %*% covBeta %*% t(X)))

MyData$SEy.pi <- sqrt(diag(X %*% covBeta %*% t(X) +
                           summary(M1)$sigma^2))

#using appropriate values
t1 <- qt(1-0.05/2,df =55-2)

#lower and upper limits for 95%
MyData$ci.ll <- MyData$P - t1 * MyData$SEy.ci
MyData$ci.ul <- MyData$P + t1 * MyData$SEy.ci

#for prediction intervals
MyData$pi.ll <- MyData$P - t1 * MyData$SEy.pi
MyData$pi.ul <- MyData$P + t1 * MyData$SEy.pi

#Figure 1.8
plot(x = HVS2$LatAbs,
     y = HVS2$OrbitalV,
     pch = 16, 
     cex = 0.8,
     xlab = "Absolute latitude",
     ylab = "Orbital volume")

with(MyData,{
polygon(c(LatAbs, rev(LatAbs)),
        c(ci.ll, rev(ci.ul)),
        col = grey(0.4), border = NULL,
        density =70)
     polygon(c(LatAbs, rev(LatAbs)),
             c(pi.ll, rev(pi.ul)),
             col = grey(0.5), border = NULL,
             density =20)
lines(LatAbs, P, lwd = 3)
})
##############################################################


##############################################################
# We can also use the predict function
PI <- predict(M1,
        newdata = MyData,
        se.fit = TRUE,
        interval = "predict")

CI <- predict(M1,
              newdata = MyData,
              se.fit = TRUE,
              interval = "confidence")
##############################################################



##############################################################
#Fit the multiple linear regression model
M1 <- lm(OrbitalV ~ LatAbs + CC + FM +
                    LatAbs : CC +
                    LatAbs : FM +
                    CC : FM,
         data = HVS2)

print(summary(M1), digits = 3, signif.stars = FALSE)
##############################################################



##############################################################
#Section 1.6 Finding the Optimal model
#Step function
step(M1)


#Fit model M2
M2 <- lm(OrbitalV ~ LatAbs + CC +
                    LatAbs : CC,
         data = HVS2)

print(summary(M2), digits = 3, signif.stars = FALSE)

#Fit model M3
M3 <- lm(OrbitalV ~ LatAbs + CC ,
         data = HVS2)
print(summary(M3), digits = 3, signif.stars = FALSE)
##############################################################



##############################################################
#Section 1.7 Degrees of freedom

X <- model.matrix(M2)
H <- X %*% solve(t(X) %*% X) %*% t(X)
sum(diag(H))

#Or use
sum(hatvalues(M2))
##############################################################



##############################################################
#Section 1.8 Model validation

#Figure 1.9
#residuals vs fitted values
E2 <- rstandard(M2)
F2 <- fitted(M2)
plot(x = F2, 
     y = E2, 
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0)
##############################################################




##############################################################
#Figure 1.10
#Cook distance plot
plot(cooks.distance(M2), 
     type = "h", 
     ylim = c(0, 1),
     ylab = "Cook distance values",
     xlab = "Observations",
     cex.lab = 1.5)
abline(h = 1, lwd = 2, lty = 2)
##############################################################


##############################################################
#Figure 1.11
#Normality
par(mfrow = c(1, 2), mar = c(5,5,2,2))
hist(E2, main = "", xlab = "Residuals", cex.lab = 1.5)
qqnorm(E2, main = "", cex.lab = 1.5)
qqline(E2)
##############################################################



##############################################################
#Figure 1.12
#Residuals vs each covariate in the model
#and not in the model
HVS2$E2 <- E2
MyVar   <- c("LatAbs", "CC", "FM",
            "Illuminance", "Temperature")
MyxyplotPolygon(HVS2, MyVar, "E2")
##############################################################



##############################################################
# Model residuals vs population
T1 <- lm(E2 ~ fPopulation, data = HVS2)
drop1(T1, test = "F")
##############################################################


##############################################################
#Figure 1.13
boxplot(E2 ~ fPopulation, 
        data = HVS2, 
        cex.lab = 1.5,
        xlab = "Population",
        ylab = "Standardized residuals")
abline(0, 0, lty = 2)
##############################################################



##############################################################
#Section 1.9 Model interpretation

#Figure 1.14 3-dimentional plot
NewData <- expand.grid(LatAbs = seq(from = 0.02, to = 65, length = 25),
                       CC     = seq(from = 1100, to = 1700, length = 25))

NewData$PredOrbitalVol <- predict(M2, NewData)
head(NewData)

pts <- data.frame(x = HVS2$LatAbs,
                  y = HVS2$CC,
                  z = HVS2$OrbitalV)

wireframe(PredOrbitalVol ~ LatAbs + CC, data = NewData,
           aspect = c(61/87, 0.4),
          light.source = c(0,0,10), 
          distance = .2,
          shade.colors.palette = function(irr, ref, height, w = .5)
          grey(w * irr + (1 - w) * (1 - (1-ref)^.4)),
          zlab = list(label = "Predicted orbital volume", rot = 90),
          shade = TRUE,
          scales = list(arrows = FALSE),
          pts = pts,
          panel.3d.wireframe =
          function(x, y, z,
                   xlim, ylim, zlim,
                   xlim.scaled, ylim.scaled, zlim.scaled,
                   pts,
                   ...) {
              xx <-
                  xlim.scaled[1] + diff(xlim.scaled) *
                      (pts$x - xlim[1]) / diff(xlim)
              yy <-
                  ylim.scaled[1] + diff(ylim.scaled) *
                      (pts$y - ylim[1]) / diff(ylim)
              zz <-
                  zlim.scaled[1] + diff(zlim.scaled) *
                      (pts$z - zlim[1]) / diff(zlim)

              panel.3dwire(x = x, y = y, z = z,
                           xlim = xlim,
                           ylim = ylim,
                           zlim = zlim,
                           xlim.scaled = xlim.scaled,
                           ylim.scaled = ylim.scaled,
                           zlim.scaled = zlim.scaled,
                           ...)

          })
##############################################################


##############################################################
#Using collinear variables
T1 <- lm(OrbitalV ~ LatAbs + CC + FM +Illuminance +
                    Temperature + fGender, 
         data = HVS2)
drop1(T1, test = "F")

T2 <- lm(OrbitalV ~ LatAbs + CC +FM, data = HVS2)
print(drop1(T2, test = "F"), digits = 3)
##############################################################


##############################################################
#Section 1.11
#Number of observations per population in case 
#we use mixed effect models
table(HVS2$fPopulation)

#Figure 1.15
boxplot(OrbitalV ~ fPopulation, 
        data = HVS2, 
        cex.lab = 1.5,
        xlab = "Population",
        ylab = "Orvital Vulume")
