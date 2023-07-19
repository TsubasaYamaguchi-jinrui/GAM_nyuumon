#    Beginner's Guide to GAM with R
#    Alain Zuur

#    www.highstat.com

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.



##############################################################
#Load the data
BL <- read.table(file = "HegerPierce.txt", 
                 header = TRUE)

#Inspect the results
str(BL)
names(BL)

##############################################################
#Housekeeping: load packages and support functions
library(lattice)
library(mgcv)

##############################################################
#Section 3.3 Linear regression

#Figure 3.2
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Depth",
     ylab ="Sources", 
     cex.lab = 1.5)


#Scale Depth
BL$DepthOriginal <- BL$Depth
BL$Depth <- BL$Depth/max(BL$Depth)

#Figure 3.3
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16,
     col = grey(0.5))

M1 <- lm(Sources ~ Depth, data = BL)
abline(M1, lwd = 5)

##############################################################
#Section 3.3 Polynimial regression
M2 <- lm(Sources ~ Depth + I(Depth^2) + I(Depth^3), data = BL)


#Figure 3.4
MD <- data.frame(Depth = seq(0.1, 1, length = 100))
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16,
     col = grey(0.5))
P2 <- predict(M2, newdata = MD)
lines(x = MD$Depth, y = P2, lwd = 5)



##############################################################
#Section 3.4 Linear spline regression

rhs <- function(x, TH) ifelse(x >= TH, x-TH,0)
M3 <- lm(Sources ~ Depth  + rhs(Depth, 0.2),
         data = BL)
P3 <- predict(M3, newdata = MD)


#Figure 3.6
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16,
     col = grey(0.5))
lines(x = MD$Depth, y = P3, lwd = 5)


#See matrix with covariates
model.matrix(M3)

#Compare models
AIC(M1, M2, M3)


#Use quantiles to decide on Knots
probs <- seq(0, 1, length = 11)
QD    <- quantile(BL$Depth, probs)
QD


#Fit the model
M4 <- lm(Sources ~ Depth  +
                   rhs(Depth, 0.159) +
                   rhs(Depth, 0.220) +
                   rhs(Depth, 0.281) +
                   rhs(Depth, 0.344) +
                   rhs(Depth, 0.410) +
                   rhs(Depth, 0.490) +
                   rhs(Depth, 0.567) +
                   rhs(Depth, 0.664) +
                   rhs(Depth, 0.787),
         data = BL)

model.matrix(M4)


#Figure 3.7
P4 <- predict(M4, newdata = MD)
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab = "Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16, 
     col = grey(0.5))
lines(x = MD$Depth, y = P4, lwd = 5)


##########################################################
#Section 3.5 Quadratic spline regression
rhs2 <- function(x, TH) ifelse(x >= TH, (x-TH)^2,0)

M5 <- lm(Sources ~ Depth  + I(Depth^2) +
                   rhs2(Depth, 0.159) +
                   rhs2(Depth, 0.220) +
                   rhs2(Depth, 0.281) +
                   rhs2(Depth, 0.344) +
                   rhs2(Depth, 0.410) +
                   rhs2(Depth, 0.490) +
                   rhs2(Depth, 0.567) +
                   rhs2(Depth, 0.664) +
                   rhs2(Depth, 0.787),
         data = BL)
        

#Figure 3.8
P5 <- predict(M5, newdata = MD)
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16,
     col = grey(0.5))
lines(x = MD$Depth, y = P5, lwd = 5)

#Numerical output
print(summary(M5), signif.stars = FALSE, digits =2)



#####################################################
#Section 3.6: Cubic regression splines

#Figure 3.9
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16, 
     col = grey(0.5))

for (i in 1:11){ abline(v = QD[i], lty = 2) }

for (i in 2:11){
 M2 <- lm(Sources ~ Depth + I(Depth^2) + I(Depth^3),
          data = BL,
          subset = BL$Depth > QD[i-1] & BL$Depth <= QD[i])

 MD <- data.frame(Depth = seq(QD[i-1], QD[i], length =100))
 P2 <- predict(M2, newdata = MD)
 lines(MD$Depth, P2, lwd = 5)
}
#End of Figure 3.9



#Fit a smoothing spline
rk<-function(x,z){
((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4 -((abs(x-z)-0.5)^4 -0.5*(abs(x-z)-0.5)^2 +7/240)/24
}


spl.X<-function(x,xk){
 q<-length(xk)+2
 n<-length(x)
 X<-matrix(1,n,q)
 X[,2]<-x
 X[,3:q]<-outer(x,xk,FUN=rk)
 X
}

X  <- spl.X(BL$Depth, QD[2:10])

print(head(X), digits=2)

M6 <- lm(Sources ~ X -1, data = BL)
summary(M6)

#Figure 3.10
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab = "Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16, 
     col = grey(0.5))

for (i in 1:11){ abline(v = QD[i], lty = 2) }

MD <- data.frame(Depth = seq(QD[1], QD[11], length =100))
Xp <- spl.X(MD$Depth, QD[2:10])
lines(MD$Depth, Xp %*% coef(M6), lwd = 5)
#End of Figure 3.10


#####################################################
#Section 3.7 The number of knots

#Fitting smoothers for a quadratic spline regression model
#using different number of knots
probs1 <- seq(0, 1, length = 3)
QD1    <- quantile(BL$Depth, probs1)

probs3 <- seq(0, 1, length = 5)
QD3    <- quantile(BL$Depth, probs3)

probs5 <- seq(0, 1, length = 5)
QD5    <- quantile(BL$Depth, probs3)

probs7 <- seq(0, 1, length = 7)
QD7    <- quantile(BL$Depth, probs7)

probs9 <- seq(0, 1, length = 9)
QD9    <- quantile(BL$Depth, probs9)

probs11 <- seq(0, 1, length = 11)
QD11    <- quantile(BL$Depth, probs11)

probs13 <- seq(0, 1, length = 12)
QD13    <- quantile(BL$Depth, probs13)

probs15 <- seq(0, 1, length = 15)
QD15    <- quantile(BL$Depth, probs15)

probs17 <- seq(0, 1, length = 17)
QD17    <- quantile(BL$Depth, probs17)

probs19 <- seq(0, 1, length = 19)
QD19    <- quantile(BL$Depth, probs19)

probs21 <- seq(0, 1, length = 21)
QD21    <- quantile(BL$Depth, probs21)


probs31 <- seq(0, 1, length = 31)
QD31    <- quantile(BL$Depth, probs31)

M5.1 <- lm(Sources ~ Depth  + I(Depth^2) +
                     rhs2(Depth, QD1[2]),
           data = BL)

M5.3 <- lm(Sources ~ Depth  + I(Depth^2) +
                     rhs2(Depth, QD3[2]) +
                     rhs2(Depth, QD3[3]) +
                     rhs2(Depth, QD3[4]),
           data = BL)

M5.5 <- lm(Sources ~ Depth  + I(Depth^2) +
                     rhs2(Depth, QD5[2]) +
                     rhs2(Depth, QD5[3]) +
                     rhs2(Depth, QD5[4]) +
                     rhs2(Depth, QD5[5]),
           data = BL)


M5.7 <- lm(Sources ~ Depth  + I(Depth^2) +
                     rhs2(Depth, QD7[2]) +
                     rhs2(Depth, QD7[3]) +
                     rhs2(Depth, QD7[4]) +
                     rhs2(Depth, QD7[5]) +
                     rhs2(Depth, QD7[6]),
           data = BL)

M5.9 <- lm(Sources ~ Depth  + I(Depth^2) +
                     rhs2(Depth, QD9[2]) +
                     rhs2(Depth, QD9[3]) +
                     rhs2(Depth, QD9[4]) +
                     rhs2(Depth, QD9[5]) +
                     rhs2(Depth, QD9[6]) +
                     rhs2(Depth, QD9[7]),
           data = BL)

M5.11 <- lm(Sources ~ Depth  + I(Depth^2) +
                      rhs2(Depth, QD11[2]) +
                      rhs2(Depth, QD11[3]) +
                      rhs2(Depth, QD11[4]) +
                      rhs2(Depth, QD11[5]) +
                      rhs2(Depth, QD11[6]) +
                      rhs2(Depth, QD11[7]) +
                      rhs2(Depth, QD11[8]),
            data = BL)

M5.13 <- lm(Sources ~ Depth  + I(Depth^2) +
                      rhs2(Depth, QD13[2]) +
                      rhs2(Depth, QD13[3]) +
                      rhs2(Depth, QD13[4]) +
                      rhs2(Depth, QD13[5]) +
                      rhs2(Depth, QD13[6]) +
                      rhs2(Depth, QD13[7]) +
                      rhs2(Depth, QD13[8]) +
                      rhs2(Depth, QD13[9]) +
                      rhs2(Depth, QD13[10]) +
                      rhs2(Depth, QD13[11]) +
                      rhs2(Depth, QD13[12]),
           data = BL)

M5.15 <- lm(Sources ~ Depth  + I(Depth^2) +
                      rhs2(Depth, QD15[2]) +
                      rhs2(Depth, QD15[3]) +
                      rhs2(Depth, QD15[4]) +
                      rhs2(Depth, QD15[5]) +
                      rhs2(Depth, QD15[6]) +
                      rhs2(Depth, QD15[7]) +
                      rhs2(Depth, QD15[8]) +
                      rhs2(Depth, QD15[9]) +
                      rhs2(Depth, QD15[10]) +
                      rhs2(Depth, QD15[11]) +
                      rhs2(Depth, QD15[12]) +
                      rhs2(Depth, QD15[13]) +
                      rhs2(Depth, QD15[14]),
            data = BL)

M5.17 <- lm(Sources ~ Depth  + I(Depth^2) +
                      rhs2(Depth, QD17[2]) +
                      rhs2(Depth, QD17[3]) +
                      rhs2(Depth, QD17[4]) +
                      rhs2(Depth, QD17[5]) +
                      rhs2(Depth, QD17[6]) +
                      rhs2(Depth, QD17[7]) +
                      rhs2(Depth, QD17[8]) +
                      rhs2(Depth, QD17[9]) +
                      rhs2(Depth, QD17[10]) +
                      rhs2(Depth, QD17[11]) +
                      rhs2(Depth, QD17[12]) +
                      rhs2(Depth, QD17[13]) +
                      rhs2(Depth, QD17[14]) +
                      rhs2(Depth, QD17[15]) +
                      rhs2(Depth, QD17[16]),
            data = BL)


M5.19 <- lm(Sources ~ Depth  + I(Depth^2) +
                      rhs2(Depth, QD19[2]) +
                      rhs2(Depth, QD19[3]) +
                      rhs2(Depth, QD19[4]) +
                      rhs2(Depth, QD19[5]) +
                      rhs2(Depth, QD19[6]) +
                      rhs2(Depth, QD19[7]) +
                      rhs2(Depth, QD19[8]) +
                      rhs2(Depth, QD19[9]) +
                      rhs2(Depth, QD19[10]) +
                      rhs2(Depth, QD19[11]) +
                      rhs2(Depth, QD19[12]) +
                      rhs2(Depth, QD19[13]) +
                      rhs2(Depth, QD19[14]) +
                      rhs2(Depth, QD19[15]) +
                      rhs2(Depth, QD19[16]) +
                      rhs2(Depth, QD19[17]) +
                      rhs2(Depth, QD19[18]),
            data = BL)

M5.31 <- lm(Sources ~ Depth  + I(Depth^2) +
                      rhs2(Depth, QD31[2]) +
                      rhs2(Depth, QD31[3]) +
                      rhs2(Depth, QD31[4]) +
                      rhs2(Depth, QD31[5]) +
                      rhs2(Depth, QD31[6]) +
                      rhs2(Depth, QD31[7]) +
                      rhs2(Depth, QD31[8]) +
                      rhs2(Depth, QD31[9]) +
                      rhs2(Depth, QD31[10]) +
                      rhs2(Depth, QD31[11]) +
                      rhs2(Depth, QD31[12]) +
                      rhs2(Depth, QD31[13]) +
                      rhs2(Depth, QD31[14]) +
                      rhs2(Depth, QD31[15]) +
                      rhs2(Depth, QD31[16]) +
                      rhs2(Depth, QD31[17]) +
                      rhs2(Depth, QD31[18]) +
                      rhs2(Depth, QD31[19]) +
                      rhs2(Depth, QD31[20]) +
                      rhs2(Depth, QD31[21]) +
                      rhs2(Depth, QD31[22]) +
                      rhs2(Depth, QD31[23]) +
                      rhs2(Depth, QD31[24]) +
                      rhs2(Depth, QD31[25]) +
                      rhs2(Depth, QD31[26]) +
                      rhs2(Depth, QD31[27]) +
                      rhs2(Depth, QD31[28]) +
                      rhs2(Depth, QD31[29]) +
                      rhs2(Depth, QD31[30]),
            data = BL)

P5.1 <- predict(M5.1, newdata = MD)
P5.3 <- predict(M5.3, newdata = MD)
P5.5 <- predict(M5.5, newdata = MD)
P5.7 <- predict(M5.7, newdata = MD)
P5.9 <- predict(M5.9, newdata = MD)
P5.11 <- predict(M5.11, newdata = MD)
P5.13 <- predict(M5.13, newdata = MD)
P5.15 <- predict(M5.15, newdata = MD)
P5.17 <- predict(M5.17, newdata = MD)
P5.19 <- predict(M5.19, newdata = MD)
P5.31 <- predict(M5.31, newdata = MD)

AllP <- c(P5.1, P5.3, P5.5, P5.7, P5.9, P5.11, P5.13, P5.15, P5.17, P5.19, P5.31)
AllD <- rep(MD$Depth, 11)
AllID <- rep(c("1 inner knot", 
               "3 inner knots", 
               "5 inner knots", 
               "7 inner knots",
               "9 inner knots",
               "11 inner knots",
               "13 inner knots",
               "15 inner knots",
               "17 inner knots",
               "19 inner knots",
               "31 inner knots"), each = length(P5.3))

AllID2 <- factor(AllID, levels =
                 c("1 inner knot",
                 "3 inner knots",
                 "5 inner knots",
                 "7 inner knots",
                 "9 inner knots",
                 "11 inner knots",
                 "13 inner knots",
                 "15 inner knots",
                 "17 inner knots",
                 "19 inner knots",
                 "31 inner knots"))

#Figure 3.11
xyplot(AllP ~ AllD | AllID2,
       col = 1, 
       type = "l", 
       lwd = 3,
       xlab = list(label = "Depth", cex = 1.5),
       ylab = list(label = "Sources", cex = 1.5),
       strip = strip.custom(bg = 'white',
       par.strip.text = list(cex = 1.2)),
       )



##################################################################
#Section 3.8 Penalized quadratic spline regression

M5 <- lm(Sources ~ Depth  + I(Depth^2) +
                   rhs2(Depth, 0.159) +
                   rhs2(Depth, 0.220) +
                   rhs2(Depth, 0.281) +
                   rhs2(Depth, 0.344) +
                   rhs2(Depth, 0.410) +
                   rhs2(Depth, 0.490) +
                   rhs2(Depth, 0.567) +
                   rhs2(Depth, 0.664) +
                   rhs2(Depth, 0.787),
         data = BL)

X <- model.matrix(M5)
coef(M5)
head(X)
K <- 9
D <- diag(rep(1, 3 + K))
D[1,1]<-D[2,2]<-D[3,3] <-0


#Figure 3.12
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16, 
     col = grey(0.5))

MyLambdas <- c(0, 1,10000)
for (lambda in MyLambdas){
  #lambda <- 0.01
  Beta.lambda <- solve(t(X) %*% X + lambda^2 * D) %*%t(X) %*% BL$Sources
  print(Beta.lambda)
  MyDepth <- seq(0.1, 1, length =100)
  XPred <- model.matrix(~ 1 + MyDepth + I(MyDepth^2)+
                        rhs2(MyDepth, 0.159) +
                        rhs2(MyDepth, 0.220) +
                        rhs2(MyDepth, 0.281) +
                        rhs2(MyDepth, 0.344) +
                        rhs2(MyDepth, 0.410) +
                        rhs2(MyDepth, 0.490) +
                        rhs2(MyDepth, 0.567) +
                        rhs2(MyDepth, 0.664) +
                        rhs2(MyDepth, 0.787) )

  Yhat <- XPred %*% as.vector(Beta.lambda)
  lines(MyDepth, Yhat, lwd = 5)
}



####################################################################
#Section 3.10 Show a Cubic smoothing spling in action

#Matrix stuff
X  <- spl.X(BL$Depth, QD[2:10])

#Wood 2006: Gives S
spl.S <- function(xk){
	q <- length(xk) + 2
	S <- matrix(0, nrow = q, ncol = q)
	S[3:q, 3:q] <- outer(xk, xk, FUN = rk)
	S
}

#Penalty matrix:
D <- spl.S(QD[2:10])


#Cross-validation
K      <- 25
OCV    <- vector(length = K)
Lambda <- seq(0, 0.001, length = K)
N      <- nrow(BL)

for (j in 1:K){
  lambda <- Lambda[j]
  EP<- vector (length= N)
  for (i in 1:N){ #Apply the GAM
       Xi    <- X[-i,]
       Betai <- solve(t(Xi) %*% Xi +lambda * D)%*%
                t(Xi)%*% BL$Sources [-i]
       Yi    <- X[i,] %*% Betai
       EP[i] <- BL$Sources[i] - Yi
  }
  OCV[j] <- sum(EP^2)/N
  }


#Figure 3.14
plot(x = Lambda, 
     y = OCV,
     type = "l",
     xlab = "lambda",
     ylab = "OCV")



#For a large data set we can use
for (j in 1:K){
  lambda <- Lambda[j]
  S      <- X %*% solve(t(X) %*% X + lambda * D) %*%t(X)
  Sii    <- diag(S)
  #Beta  <- solve(t(X) %*% X + lambda * D) %*%t(X) %*% BL$Sources
  Yfit   <- S %*% BL$Sources
  E      <- BL$Sources - Yfit
  OCV[j] <- (1/N) * sum(( E/(1-Sii))^2 )
}


#Show smoother
lambda <- 2.083333e-04
Beta   <- solve(t(X) %*% X + lambda * D) %*%t(X) %*% BL$Sources
MD     <- data.frame(Depth = seq(QD[1], QD[11], length =100))
Xp     <- spl.X(MD$Depth, QD[2:10])
fhat   <- Xp %*% Beta

#Figure 3.15
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16,
     col = grey(0.5))
lines(x = MD$Depth, y =fhat, lwd = 5)


#########################################################
#section 3.12 Degrees of freedom of a smoother

#Make a picture that visualises the relationship
#between lambda and degrees of freedom
#Use a cubic smoothing spline

X        <- spl.X(BL$Depth, QD[2:10])
D        <- spl.S(QD[2:10])   #Penalty matrix
lambda   <- 2.083333e-04
S.lambda <- X %*% solve(t(X) %*% X + lambda * D) %*%t(X)
df1      <- sum(diag(S.lambda))
df1

#df.Error <- nrow(BL) - sum(diag( 2 * S.lambda - S.lambda %*% t(S.lambda)))
#df3 <- sum(diag(S.lambda %*% t(S.lambda)))
#c(df1, df.Error, df3, nrow(BL))
###########

df     <- vector(length = 100)
Lambda <- seq(0,5, length = 100)
for (i in 1:100){
	lambda   <- Lambda[i]
    S.lambda <- X %*% solve(t(X) %*% X + lambda * D) %*%t(X)
	df[i]    <- sum(diag(S.lambda))
}

#Figure 3.16
plot(x = Lambda, 
     y = df, 
     type = "l",
     xlab = "lambda",
     ylab = "Effective degrees of freedom")


#########################################################
#Section 3.14 Confidence intervals.

#Sketch smoother with +/- 2 * SE
lambda <- 2.083333e-04
Beta   <- solve(t(X) %*% X + lambda * D) %*%t(X) %*% BL$Sources
Xp     <- spl.X(MD$Depth, QD[2:10])
yfit   <- Xp %*% Beta

S.lambda <- X %*% solve(t(X) %*% X + lambda * D) %*%t(X)
df.Error <- nrow(BL) - sum(diag( 2 * S.lambda - S.lambda %*% t(S.lambda)))
sigma2   <- sum((BL$Sources - X %*% Beta)^2) / (df.Error)
sigma    <- sqrt(sigma2)
sigma

Sp.lambda <- Xp %*% solve(t(Xp) %*% Xp + lambda * D) %*%t(Xp)
CovPredY  <- sigma2 * Sp.lambda%*% t(Sp.lambda)
SEPredY   <- sqrt(diag(CovPredY))

#Figure 3.17
par(mfrow = c(1,1), mar = c(5,5,3,3))
plot(x = BL$Depth, 
     y = BL$Sources,
     xlab = "Scaled depth",
     ylab ="Sources", 
     cex.lab = 1.5,
     cex = 0.7, 
     pch = 16,
     col = grey(0.5))
lines(x = MD$Depth, y = yfit, lwd = 3)
lines(x = MD$Depth, y = yfit + 2 * SEPredY, lwd = 5)
lines(x = MD$Depth, y = yfit - 2 * SEPredY, lwd = 5)


cbind(yfit - 2 * SEPredY,yfit, yfit + 2 * SEPredY)



###########################################################
#Section 3.15 using the function gam in mgcv
G1 <- gam(Sources ~ s(Depth), data = BL)
summary(G1)

plot(G1)

#Figure 3.18
plot(G1, 
     shift = 16.7455, 
     ylim = c(-10,100), 
     lwd =3, 
     rug = FALSE, 
     cex.lab = 1.5)
points(x = BL$Depth, 
       y = BL$Sources, 
       cex = 0.5, 
       pch = 16, 
       col = grey(0.5))


#Penalized regression regression spline
gam(Sources ~ s(Depth), data = BL)
gam(Sources ~ s(Depth,fx = TRUE,k = 5), data = BL)


G2 <-  gam(Sources ~  s(Depth, bs = "cr"), data = BL)
G2A <- gam(Sources ~ s(Depth, bs = "cr", k =20), data = BL)
G2B <- gam(Sources ~ s(Depth, bs = "cr", k =7,fx = TRUE), data = BL)
G2C <- gam(Sources ~ s(Depth, bs = "cr", k =20), gamma = 1.4, data = BL)

fitted(G1)
coefficients(G1)
X  <- predict(G1,type = "lpmatrix")
mu <- X %*% coef(G1)
X%*% G1$Vp %*% t(X)


############################################################
#Section 3.16 Danger of using Gam

E1 <- resid(G1)
F1 <- fitted(G1)

#Figure 3.19
plot(x = F1, 
     y = E1,
     cex = 0.7, 
     pch = 16,
     xlab = "Fitted values",
     ylab = "Residuals",cex.lab = 1.5)
abline(h = 0, lty = 2)


#Figure 3.20
boxplot(E1 ~ Station,  
        data = BL,
        xlab = "Station",
        ylab = "Residuals")
abline(h = 0, lty = 2)


#Figure 3.21
xyplot(E1 ~ Depth | factor(Station), 
       data = BL,
       xlab = list(label = "Depth", cex = 1.5),
       ylab = list(label = "Residuals", cex = 1.5),
         strip = strip.custom(bg = 'white',
            par.strip.text = list(cex = 1.2)),
       panel = function(x,y){
       	panel.points(x, y, pch = 16, cex = 0.7, col = 1)
       	panel.abline(h = 0, lty = 2)
       }
       )


#Figure 3.22
xyplot(Sources ~ Depth | factor(Eddy),
       groups = factor(Station),
       type = "l", 
       col = 1, 
       pch = 16,
       data = BL,
       xlab = list(label = "Depth", cex = 1.5),
       ylab = list(label = "Sources", cex = 1.5))


#Compare 2 models using AIC
G3 <- gam(Sources ~ s(Depth) + factor(Station), 
          data = BL)

G4 <- gam(Sources ~ s(Depth, by = factor(Eddy)) + factor(Station), 
          data = BL)
AIC(G3, G4)  #Small differences due to different mgcv versions


#Figure 3.23
par(mfrow = c(1,2))
plot(G4)


##########################################################
#Section 3.17 GAM with multiple smoothers

#Figure 3.24
par(mfrow = c(1,2), mar = c(5,5,3,3))
plot(x = BL$Depth, 
     y = BL$flcugl,
     xlab = "Depth",
     ylab = "flcugl",
     cex.lab = 1.5)

plot(x = BL$flcugl, 
     y = BL$Sources,
     xlab = "flcugl",
     ylab = "Sources",
     cex.lab = 1.5)


#Cubic smoothing splines
BL2 <- BL[BL$flcugl < 0.03, ]

G5 <- gam(Sources ~ s(Depth, bs = "cr") + s(flcugl, bs = "cr") + factor(Station),
          data = BL2)

#Figure 3.25
par(mfrow = c(1,2))
plot(G5)
summary(G5)


#Compare all models using AIC
G6 <- gam(Sources ~ s(Depth, bs = "cr"), data = BL2)
G7 <- gam(Sources ~ factor(Station) + 
                    s(Depth, by = factor(Eddy), bs = "cr") + 
                    s(flcugl, bs = "cr"), 
          data = BL2)
G8 <- gam(Sources ~ factor(Station) + 
                    s(Depth, by = factor(Eddy)) + 
                    s(flcugl, by = factor(Eddy), bs = "cr"), 
          data = BL2)
G9 <- gam(Sources ~ factor(Station) + 
                    s(Depth, by = factor(Eddy)) + 
                    flcugl, 
          data = BL2)
G10 <- gam(Sources ~ factor(Station) + 
                     s(Depth, by = factor(Eddy)) + 
                     flcugl * factor(Eddy), 
           data = BL2)
G11 <- gam(Sources ~ factor(Station) + te(Depth, flcugl), 
           data = BL2)
G12 <- gam(Sources ~ factor(Station) + te(Depth, flcugl, by = Eddy), 
           data = BL2)

AIC(G5, G6, G7, G8, G9, G10, G11, G12)



#Figure 3.26
par(mfrow = c(2,2))
plot(G8, scale = FALSE)

######### End of code###################
