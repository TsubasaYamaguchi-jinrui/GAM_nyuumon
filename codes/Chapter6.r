#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


##############################################################
#Load the data
Gannets <- read.table(file = "Gannets2.txt", 
                      header = TRUE)

#Inspect the results
str(Gannets)
names(Gannets)


#House keeping: load packages and support functions
library(lattice) #For multi-panel graphs
library(gamlss)  #For GAMs
library(mgcv)    #For GAMs
library(gstat)   # For checking spatial patterns

##############################################################
#Section 6.3 Brainstorming

#Figure 6.1 Spatial positions

Gannets$Xkm <- Gannets$X / 1000
Gannets$Ykm <- Gannets$Y / 1000

xyplot(Ykm ~ Xkm | factor(Year), 
       aspect = "iso",
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       data = Gannets, 
       pch=16, 
       cex = 0.4,
       col = 1,
       xlab = list(label = "Xkm", cex = 1.5),
       ylab = list(label = "Ykm", cex = 1.5))
       

#Figure 6.2
#Dotchart of Area
plot(x = Gannets$Area_surveyedkm2,
     y = 1:nrow(Gannets),
     pch = 16, 
     cex = 0.8,
     xlab = "Area",
     ylab = "Order of the data")
     
     
#Create new variables using Time
Gannets2 <- Gannets[Gannets$Area_surveyedkm2 > 0,]
Gannets2$LArea <- log(Gannets2$Area_surveyedkm2)
Gannets2$TimeH <- Gannets2$Hours + Gannets2$Minutes/60


Gannets2$Date <- paste(Gannets2$Day,Gannets2$Month,Gannets2$Year, sep = "/")
Gannets2$DayInYear <- strptime(Gannets2$Date, "%d/%m/%Y")$yday + 1

Gannets2$DaySince0 <- ceiling(
julian(strptime(Gannets2$Date, format = "%d/%m/%Y"), origin = as.Date("1991-01-01")))



##############################################################
#Section 6.4 Data exploration

#Figure 6.3
#Gannet density vs spatial coordinates and year

Gannets2$G     <- Gannets2$Gannets_in_transect
Gannets2$GProp <- Gannets2$G / Gannets2$Area_surveyedkm2

MyCex <- 2 *sqrt(20 *Gannets2$GProp / max(Gannets2$GProp))
xyplot(Y ~ X | factor(Year), 
       aspect = "iso",
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       data = Gannets2, 
       pch = 16,
       col = 1,
       layout = c(4,3),
       cex = MyCex,
       xlab = list(label = "Xkm", cex = 1.5),
       ylab = list(label = "Ykm", cex = 1.5))
       
       
#Figure 6.4
#Gannet density vs hour, day and year
xyplot(DayInYear ~ TimeH | factor(Year),
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       data = Gannets2, 
       pch = 16,
       col = 1,
       layout = c(4,3),
       cex = MyCex,
       xlab = list(label = "Hour", cex = 1.5),
       ylab = list(label = "Day", cex = 1.5))


#Figure 6.5
#Dotplot of abundance

par(mar = c(5,5,2,2))
plot(x = Gannets$G,
     y = 1:nrow(Gannets),
     xlab = "Gannets",
     ylab = "Order of the data",
     pch = 16, 
     cex = 0.8,
     cex.lab = 1.5)



##############################################################
#Section 6.5 Building up the GAM
#Fit model M1

M1 <- gam(G ~ s(Year)+ offset(LArea), 
          family = poisson,
          data = Gannets2)

E1 <- resid(M1, type = "pearson")
Overdispersion <- sum(E1^2) / M1$df.res
Overdispersion


#Fit model M2 with categorical variables
M2 <- gam(G ~ s(Year) + factor(Seastate)+ offset(LArea), 
          family = poisson,
          data = Gannets2)


#Fit Model M3 with TimeH
M3 <- gam(G ~ s(Year) + factor(Seastate) +
              s(TimeH) + offset(LArea), 
          family = poisson,
          data = Gannets2)

#Fit Model M4 with Dayin Year
M4 <- gam(G ~ s(Year) + factor(Seastate) +
              s(TimeH) + s(DayInYear) + offset(LArea), 
          family = poisson,
          data = Gannets2)

#Fit Model M5 with 2-dimentional smoothers
M5 <- gam(G ~ s(Year) + s(TimeH) + s(DayInYear)+ te(Xkm, Ykm) +
              factor(Seastate)+ offset(LArea),
         family = poisson,
         data = Gannets2)


#Fit Model M6 and put a limit in the Df
M6 <- gam(G ~ s(TimeH) + s(DayInYear) + factor(Seastate) +
              te(Xkm, Ykm, Year, k= 30) + offset(LArea), 
          family = poisson,
          data = Gannets2)


#Fit Model M7 with two 2-dimentional smoothers
M7 <- gam(G ~ s(Year) + factor(Seastate) +
              te(TimeH,DayInYear) +
              te(Xkm, Ykm) + offset(LArea), 
          family = poisson,
          data = Gannets2)



#Residuals vs Fitted values (M6)
#Figure 6.6
E6 <- resid(M6, type = "pearson")
F6 <- fitted(M6, type = "response")
plot(x = F6, 
     y = E6,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, v = 0)


#Fit model M8
M8 <- gam(G ~ s(Year) + te(TimeH,DayInYear)+
              factor(Seastate) + offset(LArea)+ te(Xkm, Ykm),
              family=negbin(c(0.1,1.5), link = log),
              data = Gannets2)
              
#Replace negbin(c(0.1,1.5)  by nb()
print(summary(M8), digits = 3, signif.stars = FALSE)


#Model validation
#Figure 6.7
par(mfrow = c(1,1), mar = c(5,5,3,3))
E8 <- resid(M8, type = "pearson")
F8 <- fitted(M8, type = "response")

MyPch <- rep(1, length = nrow(Gannets2))
MyPch[Gannets2$G > 25] <- 16
MyCex <- rep(0.5, length = nrow(Gannets2))
MyCex[Gannets2$G > 25] <- 1.2

MyCol <- rep(grey(0.4), length = nrow(Gannets2))
MyCol[Gannets2$G > 25] <- grey(0.2)


plot(x = F8, 
     y = E8,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5,
     pch = MyPch,
     cex = MyCex,
     col = MyCol)
abline(h = 0, v = 0)



#Variogram
#Figure 6.8
library(sp)
mydata <- data.frame(E8,Gannets2$Xkm, Gannets2$Ykm)
coordinates(mydata) <- c("Gannets2.Xkm", "Gannets2.Ykm")

V12 <- variogram(E8 ~ 1, mydata, cutoff = 5, robust = TRUE)
plot(x = V12$dist, 
     y = V12$gamma,
     xlab = "Distance (km)", 
     ylab = "Semi-variogram", 
     pch = 16,
     cex = 2 *V12$np / max(V12$np))



#Sketch the models
#Figure 6.9
plot(M8, select = c(1), cex.lab = 1.5)



#Figure 6.10
plot(M8, select = c(2), cex.lab = 1.5, pers = TRUE)


#Figure 6.11
plot(M8, select = c(3), cex.lab = 1.5, pers = TRUE)



##############################################################
#Section 6.6 Zero inflated GAM

#Fit the equivalent of M8
M9 <- gamlss(G ~ cs(Year, df = 8) +
                 cs(TimeH, df = 8) +
                 cs(DayInYear, df = 8) +
                 ga(~s(Xkm, Ykm, fx = TRUE, k = 28)) +
                 factor(Seastate) + offset(LArea),
              family = PO(),
              data = Gannets2)


con <- gamlss.control(n.cyc = 200)
M10 <- gamlss(G ~  cs(Year, df = 8) +
                   cs(TimeH, df = 8) +
                   cs(DayInYear,df = 8)+
                   ga(~s(Xkm, Ykm, fx = TRUE, k = 28)) +
                   factor(Seastate) + offset(LArea),
               family = ZIP(),
               data = Gannets2,
               control = con)
                   
#Extract the estimated Pi
Pi <- fitted (M10,"sigma")[1]
Pi
 
#Calculate overdispersion 
pi   <- M10$sigma.fv[1]
mu   <- M10$mu.fv
ExpY <- (1 - pi) * mu
varY <-  mu * (1 - pi) * (1 + mu * pi)
PearsonRes <- (s2$G - ExpY) / sqrt(varY)
N    <- nrow(Gannets2)
p    <- M10$df.fit
Overdispersion <- sum(PearsonRes^2) / (N - p)
Overdispersion


#Fit Model M11
M11 <- gamlss(G ~  cs(Year, df = 8) +
                   ga(~te(TimeH, DayInYear, fx = TRUE, k = 28))+
                   ga(~s(Xkm, Ykm, fx = TRUE, k = 28)) +
                   factor(Seastate) + offset(LArea),
              family = ZIP(),
              data = Gannets2)

############# End of code #######################