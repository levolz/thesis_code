############################################################################### 
## R code to produce Figure 2 for Rhemtulla, M., van Bork, R., & Borsboom, D.  
## (under review at Psych Methods)
## Written by Mijke Rhemtulla
## email: mrhemtulla@ucdavis.edu if you have comments/questions
############################################################################### 

# preliminary stuff ###########################################################
library(lavaan)
source("RhemtullaVanBorkBorsboom_Functions.R")
windowsFonts(A=windowsFont("Verdana"))
###############################################################################

# check algebra: ##############################################################

makeSigma <- function(p, b, c){
  #p = no. indicators
  #b = indicator, Y correlation
  #c = cov among indicators
  sigma <- matrix(NA, p+1, p+1)
  sigma[1:p, 1:p] <- c
  sigma[p+1, ] <- sigma[,p+1] <- b
  diag(sigma) <- 1
  
  sumMat <- matrix(0, p+1, p+2)
  sumMat[1:(p+1),1:(p+1)] <- diag(p+1)
  sumMat[1:p,(p+2)] <- 1
  
  sigma <- t(sumMat) %*% sigma %*% sumMat  
  colnames(sigma) <- c(paste("m", 1:p, sep = ""), "Y", "M")
  
  return (list(S = sigma, R = cov2cor(sigma)))
}

#what happens to cor(M,Y) as p grows? 
corvec <- NULL
for (p in 1:50){
  corvec[p] <- makeSigma(p, b = .3, c = .5)[[2]][(p+1),(p+2)]
}

plot(1:50, corvec, type = "l") #ok! 

p <- 7  #no. indicators
b <- .12 #indicator, Y correlation
c <- .48 #cov among indicators
makeSigma(p = p, b = b, c = c)

varM <- p + p*(p-1)*c #yes
varY <- 1
covMY <- p*b #yes
corMY <- covMY / (sqrt(varM)*sqrt(varY)) #yes
corMY <- p*b / sqrt(p + p*(p-1)*c) #yes

#now fit LV correlation model to sigma
sigma <- makeSigma(p, b, c)[[1]]
LVModB <- 'LM =~ m1 + m2 + m3 + m4
           LY =~ 1*Y'
fitLVB <- sem(model = LVModB, sample.cov = sigma, sample.nobs = 100000)  
inspect(fitLVB, "coef")

covMstY <- b #cov (M*, Y) 
corMstY <- b/sqrt(c)

multiplier <- sqrt(((1+(p-1)*c))/(c*p))
corMY * multiplier #YES
###############################################################################

# extend simple algebra to case w/ measurement error ##########################
makeSigma <- function(p, b, c, r){
  #p = no. indicators
  #b = indicator, Y correlation
  #c = cov among latent causal indicators (m)
  #r = reliability of causal indicators 
  sigma <- matrix(NA, 2*p+1, 2*p+1)
  sigma[1:p, 1:p] <- c
  sigma[(p+1):(2*p),(p+1):(2*p)] <- r*c
  sigma[(p+1):(2*p),1:p] <- sigma[1:p,(p+1):(2*p)] <- sqrt(r)*c
  if (p == 1){
    sigma[1,2] <- sigma[2,1] <- sqrt(r)
  } else {
    diag(sigma[(p+1):(2*p),1:p]) <- diag(sigma[1:p,(p+1):(2*p)]) <- sqrt(r)
  }
  sigma[(2*p)+1,1:p] <- sigma[1:p,(2*p)+1] <- b
  sigma[(2*p)+1,(p+1):(2*p)] <- sigma[(p+1):(2*p),(2*p)+1] <- sqrt(r)*b
  diag(sigma) <- 1
  
  sumMat <- matrix(0, (2*p)+1, (2*p)+3)
  sumMat[1:(2*p+1),1:(2*p+1)] <- diag((2*p)+1)
  sumMat[1:p,(2*p+2)] <- 1
  sumMat[(p+1):(2*p),(2*p+3)] <- 1
  
  sigma <- t(sumMat) %*% sigma %*% sumMat  
  colnames(sigma) <- c(paste("m", 1:p, sep = ""), paste("x", 1:p, sep = ""), "Y", "M", "X")
  
  return (list(S = sigma, R = cov2cor(sigma)))
}

#what happens to cor(M,Y) as p grows? 
corMYvec <- NULL
corXYvec <- NULL
for (p in 1:50){
  corMYvec[p] <- makeSigma(p, b = .3, c = .5, r = .8)[[2]][(2*p+1),(2*p+2)]
  corXYvec[p] <- makeSigma(p, b = .3, c = .5, r = .8)[[2]][(2*p+1),(2*p+3)]
}

plot(1:50, corMYvec, type = "l") #ok! 
lines(1:50, corXYvec, col = "red") #ok! 

p <- 7  #no. indicators
b <- .12 #indicator, Y correlation
c <- .48 #cov among latent indicators
r <- .7  #reliability of observed indicators
makeSigma(p = p, b = b, c = c, r = r)

varM <- p + p*(p-1)*c #yes
varX <- p + p*(p-1)*c*r #yes
varY <- 1

covMY <- p*b #yes
corMY <- covMY / (sqrt(varM)*sqrt(varY)) #yes
corMY <- p*b / sqrt(p + p*(p-1)*c) #yes

covMX <- p*sqrt(r) + p*(p-1)*sqrt(r)*c #yes!
corMX <- covMX / (sqrt(varX)*sqrt(varM)) #yes

covXY <- sqrt(r)*b*p #yes
corXY <- (sqrt(r)*b*p) / sqrt(p + p*(p-1)*c*r) #yes


#now fit LV correlation model to sigma
sigma <- makeSigma(p, b, c, r)[[1]]
LVModB <- 'LM =~ x1 + x2 + x3 + x4
LY =~ 1*Y'
fitLVB <- sem(model = LVModB, sample.cov = sigma, sample.nobs = 100000)  
inspect(fitLVB, "coef")

varMst <- r*c #yes
covMstY <- sqrt(r)*b #cov (M*, Y)  #yes b/c this is the total cov(x_i, Y) and lambda = 1
inspect(fitLVB, "coef")$psi[1,2] #yes this is sqrt(r)*b

#std model: 
lambda <- sqrt(r*c) #yes
corMstY <- sqrt(r)*b/sqrt(r*c) #yes!
corMstY <- b/sqrt(c) #yes!
inspect(fitLVB, "std")$psi[1,2] #yes

#figure out multiplier:  
multiplierMMst <- sqrt(((1+(p-1)*c))/(c*p))
corMY * multiplierMMst == corMstY#YES

multiplierMX <-  sqrt((r + (p-1)*c*r) / (1 + (p-1)*c*r))
corMY * multiplierMX #YES == corXY
###############################################################################

# Make Figure 2 ###############################################################
makeMult <- function(p, c, r) {
  return(list(multiplierMMst = sqrt(((1+(p-1)*c))/(c*p)) * (r/r), 
              multiplierMX = sqrt((r + (p-1)*c*r) / (1 + (p-1)*c*r))))
}

#vary item covariance 
c <- seq(0, 1, by = .01)
m3 <- makeMult(p = 3, c = c, r = .7)[[1]]
m5 <- makeMult(p = 5, c = c, r = .7)[[1]]
m8 <- makeMult(p = 8, c = c, r = .7)[[1]]
m12 <- makeMult(p = 12, c = c, r = .7)[[1]]

x3 <- makeMult(p = 3, c = c, r = .7)[[2]]
x5 <- makeMult(p = 5, c = c, r = .7)[[2]]
x8 <- makeMult(p = 8, c = c, r = .7)[[2]]
x12 <- makeMult(p = 12, c = c, r = .7)[[2]]

#vary reliability: 
r <- seq(0, 1, by = .01)
m3r <- makeMult(p = 3, c = .5, r = r)[[1]]
m5r <- makeMult(p = 5, c = .5, r = r)[[1]]
m8r <- makeMult(p = 8, c = .5, r = r)[[1]]
m12r <- makeMult(p = 12, c = .5, r = r)[[1]]

x3r <- makeMult(p = 3, c = .5, r = r)[[2]]
x5r <- makeMult(p = 5, c = .5, r = r)[[2]]
x8r <- makeMult(p = 8, c = .5, r = r)[[2]]
x12r <- makeMult(p = 12, c = .5, r = r)[[2]]

jpeg("multiplier_measErr_cr_4figs.jpg", width = 8, height = 8, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")

par(fig = c(0, .5, .5, 1))
plot (c, m3, type = "l", col = 1, lty = 1, ylim = c(1, 2),
      xlab = "correlations among indicators", 
      ylab = "LV bias")  
lines(c, m5, type = "l", col = 2, lty = 2)  
lines(c, m8, type = "l", col = 3, lty = 3)  
lines(c, m12, type = "l", col = 4, lty = 4)  
legend(x = .5, y = 2, legend = c("p = 3", "p = 5", "p = 8", "p = 12"),
       lty = 1:4, col = 1:4, bty = "n")

par(fig = c(.5, 1, .5, 1), new = TRUE)
plot (c, x3, type = "l", col = 1, lty = 1, ylim = c(.5, 1),
      xlab = "correlations among indicators", 
      ylab = "composite bias")  
lines(c, x5, type = "l", col = 2, lty = 2)  
lines(c, x8, type = "l", col = 3, lty = 3)  
lines(c, x12, type = "l", col = 4, lty = 4)  
legend(x = .5, y = .7, legend = c("p = 3", "p = 5", "p = 8", "p = 12"),
       lty = 1:4, col = 1:4, bty = "n")

par(fig = c(0, .5, 0, .5), new = TRUE)
plot(r, m3r, type = "l", col = 1, lty = 1,
     xlab = "item reliability", ylab = "LV bias",
     xlim = c(0, 1), ylim = c(1, 2))  
lines(r, m5r, type = "l", col = 2, lty = 2)  
lines(r, m8r, type = "l", col = 3, lty = 3)  
lines(r, m12r, type = "l", col = 4, lty = 4)  
legend(x = .5, y = 2, legend = c("p = 3", "p = 5", "p = 8", "p = 12"),
       lty = 1:4, col = 1:4, bty = "n")

par(fig = c(.5, 1, 0, .5), new = TRUE)
plot(r, x3r, type = "l", col = 1, lty = 1,
     xlab = "item reliability", ylab = "composite bias",
     xlim = c(0, 1), ylim = c(.5, 1))  
lines(r, x5r, type = "l", col = 2, lty = 2)  
lines(r, x8r, type = "l", col = 3, lty = 3)  
lines(r, x12r, type = "l", col = 4, lty = 4)  
legend(x = .5, y = .7, legend = c("p = 3", "p = 5", "p = 8", "p = 12"),
       lty = 1:4, col = 1:4, bty = "n")

dev.off()

###############################################################################
