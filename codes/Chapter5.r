#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


##############################################################
#Load the data
Squid <- read.table(file = "SquidNorway.txt",
                    header = TRUE)

#Inspect the results
str(Squid)
names(Squid)


##############################################################
#Housekeeping: load packages and support functions
library(lattice) #For multi-panel graphs
library(mgcv)    #For GAMs

#Source our library file Highstatlib.R
source(file = "HighstatLibV8.R")
#It can be dowloaded from the course website.
##############################################################



##############################################################
#Section 5.3 Data exploration

#Figure 5.1
MyVar <- c("Lat", "Depth", "ML", "d15N")
Mydotplot(as.matrix(Squid[,MyVar]))


#Collinerity
#Figure 5.2
MyVar <- c("d15N" , "Lat", "Depth", "ML")
pairs(Squid[,MyVar], lower.panel = panel.cor)


#Multipanel scatterplot
#Figure 5.3
MyVar <- c("Lat", "Depth", "ML")
Myxyplot(Squid, MyVar, "d15N")

#For the book chapter we changed the layout
#For this we need the xyplot code, and not the
#Myxyplot function

Y   <- Squid$d15N
MyX <- c("Lat", "Depth", "ML")
X   <- Squid[, MyX]

AllX  <- as.vector(as.matrix(Squid[,MyX]))
AllY  <- rep(Squid[,"d15N"] , length(MyX))
AllID <- rep(MyX, each = nrow(Squid))
   
xyplot(AllY ~ AllX|factor(AllID), 
       col = 1,
       layout = c(3,1),
       xlab = list(label="Explanatory variables",cex = 1.5),
       ylab = list(label = "d15N", cex = 1.5),
       strip = function(bg='white', ...)
       strip.default(bg='white', ...),
       scales = list(alternating = T,
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x, y){
                panel.grid(h=-1, v= 2)
                panel.points(x, y, col = 1)
                panel.loess(x, y, col = 1, lwd = 2)})
#Using the Myxyplot is more friendly on the eye!
#Have a look at ggplot2...yiu can make similar graphs,
#but with easier code.


##############################################################
#Section 5.5 Applying the multiple linear regression model

#Fit Model M1
M1 <- glm(d15N ~ ML + Lat, data = Squid)

#Model validation
#Figure 5.4
E1 <- rstandard(M1)
F1 <- fitted(M1)

par(mfrow = c(1,2), mar = c(5,5,3,3))
plot(x = F1, 
     y = E1, 
     cex.lab = 1.5,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0)

plot(x = Squid$ML, 
     y = E1, 
     cex.lab = 1.5,
     xlab = "Mantel length",
     ylab = "Residuals")
abline(h = 0)


#Fit a GAM using E1
T1 <- gam(E1 ~ -1 + s(ML), data = Squid)
summary(T1)

#Figure 5.5
par(mar = c(5,5,3,3))
plot(T1, cex.lab = 1.5)
abline(h = 0, lty = 2)


##############################################################
#Section 5.6 Applying an additive model

#Fit model M2
M2 <- gam(d15N ~ Lat + s(ML), data = Squid)
summary(M2)

#Figure 5.6
plot(M2, resid = TRUE, pch = 16, cex = 0.5)
abline(h = 0, lty = 2)


#############################################################
#Model validation
#Figure 5.7
E2 <- resid(M2)
F2 <- fitted(M2)
par(mfrow = c(2,3), mar = c(5,5,3,3))
plot(x = F2, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2)

plot(x = Squid$ML, 
     y = E2,
     xlab = "ML",
     ylab = "Residuals")
abline(h = 0, lty = 2)

plot(x = Squid$Lat, 
     y = E2,
     xlab = "Latitude",
     ylab = "Residuals")
abline(h = 0, lty = 2)

plot(x = Squid$Depth, 
     y = E2,
     xlab = "Depth",
     ylab = "Residuals")
abline(h = 0, lty = 2)

hist(E2, xlab = "", ylab ="", breaks = 10)

plot(sort(E2), 
     type = "h",
     xlab = "Sorted residuals",
     ylab = "Residuals")
     
     
     
##############################################################
#Section 5.7 Testing linearity vs non-linearity


#Compare models
AIC(M1, M2)
anova(M1, M2, test = "F")

#Fit model M3
M3 <- gam(d15N ~ Lat + s(ML, bs = "cr"), data = Squid)


#5.7.1 Programming a smoother manually
rg <- range(Squid$ML)
Squid$MLsc <- (Squid$ML - rg[1]) / (rg[2] - rg[1])

probs <- seq(0, 1, length = 5)
QD    <- quantile(unique(Squid$MLsc), probs)
QD

rhs <- function(x, TH) {ifelse(x >= TH, (x-TH)^3,0)}
dk  <- function(x,TH,K){ (rhs(x,TH) - rhs(x,K)) / (K-TH) }
bj  <- function(x,TH,K){ dk(x,TH,K) - dk(x,K-1,K)}

I1     <- order(Squid$MLsc)
Squid1 <- Squid[I1,]

M4 <- lm(d15N ~ 1 + Lat + MLsc  +
                bj(MLsc, QD[2], QD[4]) +
                bj(MLsc, QD[3], QD[4]) ,
         data = Squid1)

X <- model.matrix(M4)
head(X)

coef(M4)
Smooth <- X[,3:5] %*% coef(M4)[3:5]
Smooth <-Smooth - mean(Smooth)


#Figure 5.8
plot(x = Squid1$MLsc, 
     y = Squid1$d15N,
     xlab = "Scaled ML",
     type = "n",
     ylab = "Smoothing function ML",
     ylim = c(-2.5,2.5),
     cex = 0.7, 
     pch = 16, 
     col = grey(0.5))

E4 <- resid(M4)
lines(Squid1$MLsc, Smooth, lwd = 5)
for (i in 1:5){abline(v = QD[i])}
points(Squid1$MLsc, Smooth+E4)
#End Figure 5.8 code



#Compare the model with lm results
M5 <- lm(d15N ~ Lat + MLsc , data = Squid1)
anova(M4, M5, test = "F")
#Small differences with the book!


#############################################################
#Section 5.8 Consequences of ignoring collinearity in GAM

M6 <- gam(d15N ~ s(Lat)+ s(Depth) + s(ML), data = Squid)
M6 <- gam(d15N ~ s(Lat, k = 4) + s(Depth, k = 4) + s(ML), data = Squid)
summary(M6)

##################End of code ##########################

