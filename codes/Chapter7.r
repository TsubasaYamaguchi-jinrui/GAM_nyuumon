#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


##############################################################
#Load the data
Par <- read.table("Parasites3.txt",
                  header = TRUE)

Par$Worms <- Par$Elytrophalloides_oatesi

#Inspect the data
str(Par)
names(Par)

##############################################################
#House keeping: load packages and support functions
library(lattice) #For multi-panel graphs
library(mgcv)    # For GAMs

#Source our library file Highstatlib.R
#source(file="HighstatLibV4.R")
#You need to adjust the path!
#It can be dowloaded from the course website.


##############################################################
#Section 7.3 Data exploration

#Sampling effort over time and space
with(Par, table(Year, Area))

#Figure 7.1
boxplot(Length ~ Sex, 
        data = Par,
        xlab = "Sex",
        ylab = "Length")


#Figure 7.2
ylab.name = expression(paste("Prevalence of ",italic("E. oatesi"), sep=""))

xyplot(Worms ~ Length | Sex * Area,
       pch = 16,
       xlab = list(label = "Length", cex = 1.5),
       ylab = list(label = ylab.name, cex = 1.5),
       data = Par, 
       col=1)
       

#table command
table(Par$Worms)
   
   
       
##############################################################
#Section 7.5 Applying Binomial GAM

M1 <- gam(Worms ~ s(Length), 
          data = Par,
          family = binomial)
         
summary(M1)

deviance(M1)
logLik(M1)
AIC(M1)


#Figure 7.3
E1 <- resid(M1, type = "pearson")
F1 <- fitted(M1, type = "response")

par(mfrow = c(1,1), mar = c(5,5,3,3))
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values",
     ylab = "Pearson residuals", 
     pch = 16)
abline(h = 0)


#Figure 7.4
plot(M1, 
     resid = TRUE, 
     shade = TRUE,
     cex = 0.5,
     pch = 16)



#Figure 7.5
range(Par$Length)
ylab.name = expression(paste("Probability of presence of ",
                             italic("E. oatesi"), 
                             sep=""))

MD <- data.frame(Length = seq(29.5, 57, length = 100))
P1 <- predict(M1, newdata = MD, se = TRUE, type = "link")

SE.UP  <- exp(P1$fit + 2 * P1$se.fit) / (1+exp(P1$fit + 2 * P1$se.fit))
SE.Low <- exp(P1$fit - 2 * P1$se.fit) / (1+exp(P1$fit - 2 * P1$se.fit))
Fit    <- exp(P1$fit) / (1+exp(P1$fit))

plot(x = Par$Length,  
     y = Par$Worms,
     xlab = "Length",
     ylab = ylab.name,
     pch = 16, col = grey(0.3))

polygon(c(MD$Length, rev(MD$Length)),
        c(SE.Low, rev(SE.UP)),
        col = grey(0.5),border=NULL,
        density = 50   )

lines(MD$Length, Fit, lwd = 3)




#Fit Model M2:
M2 <- gam(Worms ~ s(Length) + Sex, 
          data = Par,
         family = binomial)
anova(M2)
AIC(M1, M2)


#Fit Model M3
M3 <- gam(Worms ~ s(Length, by = Sex) + Sex, 
          data = Par,
          family = binomial)         
summary(M3)

AIC(M1, M2, M3)


#Figure 7.6
par(mfrow = c(1,2), mar = c(5,5,3,3))
plot(M3, 
     resid = TRUE, 
     shade = TRUE,
     cex = 0.5, 
     cex.lab = 1.5, 
     pch = 16)



#Fit model M4
M4 <- gam(Worms ~ s(Length, by = Sex) + Sex + Area, 
          data = Par,
          family = binomial)
          
AIC(M1, M2, M3, M4)
summary(M4)


#Figure 7.7
par(mfrow = c(1,2), mar = c(5,5,3,3))
plot(M4, 
     resid = TRUE, 
     shade = TRUE,
     cex = 0.5, 
     cex.lab = 1.5, 
     pch = 16)



#Figure 7.8
Raise <- (1 - Par$Worms) * as.numeric(Par$Area)/50
Lower <- Par$Worms * as.numeric(Par$Area)/50
MyPch <- rep(1, nrow(Par))
MyPch[Par$Area == "SanMatias"] <- 16


par(mfrow = c(1,1), mar = c(5,5,3,3))
plot(x = Par$Length, 
     y = Par$Worms + Raise - Lower,
     xlab = "Length",
     ylab = ylab.name,
     col = 1, 
     pch = MyPch,
     cex.lab = 1.5)

polygon(c(MD$Length, rev(MD$Length)),
        c(SE.Low, rev(SE.UP)),
        col = grey(0.5),border=NULL,
        density = 50   )

lines(x = MD$Length, y = Fit, lwd = 3)



#Fit Models M5, M6, M7
M5 <- gam(Worms ~ s(Length) + Sex + Area, 
          data = Par,
          family = binomial)

M6 <- gam(Worms ~ Length * Sex + Area, 
          data = Par,
          family = binomial)

M7 <- glm(Worms ~ Length * Sex + Area, 
          data = Par,
          family = binomial)

AIC(M1, M2, M3, M4, M5, M6, M7)
drop1(M7, test= "Chi")


#Figure 7.9
MD <- expand.grid(Length = seq(29.5, 57, length = 100),
                  Sex    = levels(Par$Sex),
                  Area   = levels(Par$Area))

P7 <- predict(M7, newdata = MD, se = TRUE, type = "link")
MD$fit <- P7$fit
MD$se.fit <- P7$se.fit


xyplot(Worms ~ Length | Area * Sex,
     data = Par,
     xlab = list(label = "Length", cex = 1.5),
     ylab = list(ylab.name, cex = 1.5),
     panel = function(x,y,subscripts,...){
     	panel.points(x, y, cex = 1, pch = 16, col = 1)
     	Sexi   <- as.character(Par[subscripts,"Sex"][1])
     	Areai  <- as.character(Par[subscripts,"Area"][1])
     	MDi    <- MD[MD$Sex==Sexi & MD$Area == Areai,]
        Fit    <- exp(MDi$fit) / (1+exp(MDi$fit))
        SE.UP  <- exp(MDi$fit + 2 * MDi$se.fit) / (1+exp(MDi$fit + 2 * MDi$se.fit))
        SE.Low <- exp(MDi$fit - 2 * MDi$se.fit) / (1+exp(MDi$fit - 2 * MDi$se.fit))
        panel.polygon(c(MDi$Length, rev(MDi$Length)),
                c(SE.Low, rev(SE.UP)),
                col = grey(0.5),border=NULL,
                density = 50   )
     	panel.lines(MDi$Length, Fit, lwd = 3, col = 1)
     		})

########## End of code##############



