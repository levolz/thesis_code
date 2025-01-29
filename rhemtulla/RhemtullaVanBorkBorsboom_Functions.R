################################################################################################### 
## R code to produce Figures for Rhemtulla, M., van Bork, R., & Borsboom, D.  
## (under review at Psych Methods)
## Written by Mijke Rhemtulla
## email: mrhemtulla@ucdavis.edu if you have comments/questions
################################################################################################### 

# 1. load libraries and lavaan model syntax ################################### 
library(lavaan)
library(R.utils)

CompMod <- 'M2 ~ a*X
            Y ~ b*M2
            ab := a*b'

LVMod <- 'LM =~ m1 + m2 + m3 + m4
          LM ~ a*X
          Y ~ b*LM
          ab := a*b'

FullMod <- 'Y ~ m1 + m2 + m3 + m4
            m1 ~ X
            m2 ~ X
            m3 ~ X
            m4 ~ X'

netMod <- 'Y ~ m4
           m4 ~ m3 + m2
           m3 ~ m1
           m2 ~ m1
           m1 ~ X'

MimicMod <- 'Y ~ LM
             LM =~ m3 + m4
             LM ~ m1 + m2
             m1 ~~ m2
             m1 ~ X
             m2 ~ X'

################################################################################################### 

#2. function to generate data for reliable compositea/causal indicator model ######################
generateCov <- function(a = NA, b = .4, avec = c(.5, .3, .2, 0), covm = NA, res.M = 0){
  #a: total influence of X on M (composite)
  #b: regression of Y on M
  #a1 - a4 are regression paths of causal indicators of M on X:
  #res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m4)
  #nb this function returns only perfectly reliable components!
  
  amat <- avec %*% t(avec)
  
  gamma <- a/sum(avec)
  covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)   #NEW so total var(M) includes residual variance
  
  #variance of residuals for the causal indicators of M:
  evec <- 1 - avec^2
        
  #residual covariances among m1-m4: 
  emat <- covm - amat  
  
  Beta <- matrix(0, 7, 7)
  Beta[1:4,7] <- avec
  Beta[5,1:4] <- gamma
  Beta[6,5] <- b
  I <- diag(7)
  Psi <- matrix(0, 7, 7)
  Psi[1:4, 1:4] <- emat
  diag(Psi) <- c(evec, res.M, (1-b^2), 1)  
   
  Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
 
  #add a column for M (sum of m1-m4)
  sumMat <- matrix(0, 7, 8)
  sumMat[1:7,1:7] <- diag(7)
  sumMat[1:4,8] <- 1
  
  Sigma <- t(sumMat) %*% Sigma %*% sumMat  
  
  colnames(Sigma) <- c("m1", "m2", "m3", "m4", "M", "Y", "X", "M2")
  return(Sigma)  
}

# TO USE
# Sigma <- generateCov(a = .4, b = .4, avec = c(.4, .4, .4, .4), covm = NA, res.M = 0)
#
# fitLV <- sem(model = LVMod, sample.cov = Sigma, sample.nobs = 100000, std.lv = TRUE)  
# fitFull <- sem(model = FullMod, sample.cov = Sigma, sample.nobs = 100000, sample.cov.rescale = FALSE)
# 
# summary(fitFull)
# summary(fitComp)
# summary(fitLV)
################################################################################################### 

#3. fx to generate data for unreliable composite/causal indicators ################################ 
#function to generate from composite/formative model with unreliable formative indicators: 
generateCov3 <- function(a = .6, b = .4, avec = c(.5, .3, .2, 0), covm = NA, res.M = 0, r = .8){
  #a: total influence of X on M (composite)
  #b: regression of Y on M
  #a1 - a4 are regression paths of causal indicators of M on X:
  #res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m4)
  #r is the reliability of the composite/causal indicators
  amat <- avec %*% t(avec)
  
  # solve for gamma and covm given a
  gamma <- a/sum(avec)
  covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)   #NEW so total var(M) includes residual variance
  
  #variance of residuals for the causal indicators of M:
  evec <- 1 - avec^2
  
  #residual covariances among m1-m4: 
  emat <- covm - amat  
  
  Beta <- matrix(0, 7, 7)
  Beta[1:4,7] <- avec
  Beta[5,1:4] <- gamma
  Beta[6,5] <- b
  I <- diag(7)
  Psi <- matrix(0, 7, 7)
  Psi[1:4, 1:4] <- emat
  diag(Psi) <- c(evec, res.M, (1-b^2), 1)  
  
  Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
  
  #create unreliable indicators of m1-m4, call them y1-y4:
  Lambda <- matrix(0, 11, 7)
  diag(Lambda[1:7, 1:7]) <- 1
  diag(Lambda[8:11, 1:4]) <- sqrt(r)
  Theta <- matrix(0, 11, 11)
  diag(Theta[8:11, 8:11]) <- 1-r
  
  Sigma2 <- Lambda %*% Sigma %*% t(Lambda) + Theta
  
  #add a column for M (sum of y1-y4)
  sumMat <- matrix(0, 11, 12)
  sumMat[1:11,1:11] <- diag(11)
  sumMat[8:11,12] <- 1
  
  Sigma3 <- t(sumMat) %*% Sigma2 %*% sumMat  
  
  colnames(Sigma3) <- c("m1", "m2", "m3", "m4", "M", "Y", "X", "y1", "y2", "y3", "y4", "M2")
  
  return(Sigma3)  
}

################################################################################################### 

#4. variation w/ unequal gammas (causal indicator loadings)######################################## 
#function to generate from composite/formative model with unequal gammas (loadings):
#NB this means the composite model (which assumes equal weights) is no longer accurate

generateCov2 <- function(a = NA, b = .4, avec = c(.5, .3, .2, .4), covm = NA, res.M = 0){
  #a: total influence of X on M (composite)
  #b: regression of Y on M
  #a1 - a4 are regression paths of causal indicators of M on X:
  #res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m4)
  amat <- avec %*% t(avec)
  
  gammamean <- a/sum(avec)
  sumsqgamma <- 1
  
  ptm <- proc.time()
  while(sumsqgamma >= 1 & (proc.time() - ptm)[3] < 30){ #quit after 30 seconds
    gamma123 <- rnorm(3, mean = gammamean, sd = .1)
    gamma4 <- (a - sum(gamma123*avec[1:3]))/avec[4]
    gammavec <- c(gamma123, gamma4)
    sumsqgamma <- sum(gammavec^2)
  }
  
  gammamat <- gammavec %*% t(gammavec)
  diag(gammamat) <- 0
  covm <- (1-sumsqgamma - res.M) / sum(gammamat)   #NEW so total var(M) includes residual variance
  
  #variance of residuals for the causal indicators of M:
  evec <- 1 - avec^2
  
  #residual covariances among m1-m4: 
  emat <- covm - amat  
  
  Beta <- matrix(0, 7, 7)
  Beta[1:4,7] <- avec
  Beta[5,1:4] <- gammavec
  Beta[6,5] <- b
  I <- diag(7)
  Psi <- matrix(0, 7, 7)
  Psi[1:4, 1:4] <- emat
  diag(Psi) <- c(evec, res.M, (1-b^2), 1)  
  
  Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
  
  #add a column for M (sum of m1-m4)
  sumMat <- matrix(0, 7, 8)
  sumMat[1:7,1:7] <- diag(7)
  sumMat[1:4,8] <- 1
  
  Sigma <- t(sumMat) %*% Sigma %*% sumMat  
  
  colnames(Sigma) <- c("m1", "m2", "m3", "m4", "M", "Y", "X", "M2")
  return(list(Sigma, gammavec, covm))  
}

################################################################################################### 

#5. fxs to generate data for MIMIC model ########################################################### 
generateCov.mimic <- function(a = NA, b = .4, avec = c(.5, .3), lamvec = c(.6, .7), covm = NA, res.M = 0){
  #a: total influence of X on M
  #b: regression of Y on M
  #a1 - a2 are regression paths of causal indicators of M on X:
  #a3 - a4 are factor loadings of reflective indicators of M
  #res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m2)
  amat <- avec %*% t(avec)
  
  #if covm is not provided, then solve for gamma and covm
  if (is.na(covm) == TRUE) {
    gamma <- a/sum(avec)
    covm <- (1 - 2*gamma^2 - res.M)/(2*gamma^2)   #NEW so total var(M) includes residual variance
  } else {
    #if a is not provided, then solve for gamma and a
    gamma <- 1/(2*covm + 2)
    a <- gamma*(sum(avec)) 
  }
  
  #variance of residuals for the causal indicators of M:
  evec <- 1 - avec^2
  
  #residual covariances among m1-m4: 
  emat <- covm - amat  
  
  Beta <- matrix(0, 7, 7)
  Beta[1:2,7] <- avec
  Beta[5,1:2] <- gamma
  Beta[3:4,5] <- lamvec
  Beta[6,5] <- b
  I <- diag(7)
  Psi <- matrix(0, 7, 7)
  Psi[1:2, 1:2] <- emat
  diag(Psi) <- c(evec, 1 - lamvec^2, res.M, (1-b^2), 1)  
  
  Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
  
  #add a column for M (sum of m1-m4)
  sumMat <- matrix(0, 7, 8)
  sumMat[1:7,1:7] <- diag(7)
  sumMat[1:4,8] <- 1
  
  Sigma <- t(sumMat) %*% Sigma %*% sumMat  
  
  colnames(Sigma) <- c("m1", "m2", "m3", "m4", "M", "Y", "X", "M2")
  return(Sigma)  
}

generateCov.mimic.r <- function(a = NA, b = .6, avec = c(.5, .3), lamvec = c(.6, .7), covm = NA, res.M = 0, r = .8){
  #a: total influence of X on M
  #b: regression of Y on M
  #a1 - a2 are regression paths of causal indicators of M on X:
  #a3 - a4 are factor loadings of reflective indicators of M
  #res.M is the residual variance of M (i.e., if it is more than just the sum of m1-m2)
  #r is the reliability of each observed causal indicator of M
  amat <- avec %*% t(avec)
  
  #if covm is not provided, then solve for gamma and covm
  if (is.na(covm) == TRUE) {
    gamma <- a/sum(avec)
    covm <- (1 - 2*gamma^2 - res.M)/(2*gamma^2)
  } else {
    #if a is not provided, then solve for gamma and a
    gamma <- 1/(2*covm + 2)
    a <- gamma*(sum(avec)) 
  }
  
  #variance of residuals for the causal indicators of M:
  evec <- 1 - avec^2
  
  #residual covariances among m1-m4: 
  emat <- covm - amat  
  
  Beta <- matrix(0, 7, 7)
  Beta[1:2,7] <- avec
  Beta[5,1:2] <- gamma
  Beta[3:4,5] <- lamvec
  Beta[6,5] <- b
  I <- diag(7)
  Psi <- matrix(0, 7, 7)
  Psi[1:2, 1:2] <- emat
  diag(Psi) <- c(evec, 1 - lamvec^2, res.M, (1-b^2), 1)  
  
  Sigma <- solve(I - Beta) %*% Psi %*% t(solve(I - Beta))
  
  #create unreliable indicators of m1-m2, call them y1-y2:
  Lambda <- matrix(0, 9, 7)
  diag(Lambda[1:7, 1:7]) <- 1
  diag(Lambda[8:9, 1:2]) <- sqrt(r)
  Theta <- matrix(0, 9, 9)
  diag(Theta[8:9, 8:9]) <- 1-r
  
  Sigma2 <- Lambda %*% Sigma %*% t(Lambda) + Theta
  
  #add a column for M (sum of y1-y2 + m3-m4)
  sumMat <- matrix(0, 9, 10)
  sumMat[1:9,1:9] <- diag(9)
  sumMat[c(3:4,8:9),10] <- 1
  
  Sigma3 <- t(sumMat) %*% Sigma2 %*% sumMat  
  
  colnames(Sigma3) <- c("m1", "m2", "m3", "m4", "M", "Y", "X", "y1", "y2", "M2")
  
  return(Sigma3)  
}
################################################################################################### 

#6. fxs to generate data for network model ######################################################## 
generateCov.net <- function(a = .4, b = .4, beta = c(.3, .3, .3, .3)){
  #total effect of X -> m4 = a: 
  a. <- a/(beta[1]*beta[3] + beta[2]*beta[4])
  
  B <- matrix(0, 6, 6)
  B[2,1] <- a.
  B[3,2] <- beta[1]
  B[4,2] <- beta[2]
  B[5,3] <- beta[3]
  B[5,4] <- beta[4]
  B[6,5] <- b
  
  I <- diag(6)
  evec <- diag(B %*% t(B))
  evec[5] <- evec[5] + 2*beta[1]*beta[2]*beta[3]*beta[4]
  
  Psi <- diag(6) - diag(evec)
  Sigma <- solve(I - B) %*% Psi %*% t(solve(I - B))
  
  #add a column for M (sum of m1-m4)
  sumMat <- matrix(0, 6, 7)
  sumMat[1:6,1:6] <- diag(6)
  sumMat[2:5,7] <- 1
  
  Sigma <- t(sumMat) %*% Sigma %*% sumMat
  
  colnames(Sigma) <- c("X", "m1", "m2", "m3", "m4", "Y", "M2")
  return(Sigma)    
}

generateCov.net.r <- function(a = .4, b = .4, beta = c(.3, .3, .3, .3), r = .8){
  #total effect of X -> m4 = a: 
  #r = reliability of observed indicators 
  a. <- a/(beta[1]*beta[3] + beta[2]*beta[4])
  
  B <- matrix(0, 6, 6)
  B[2,1] <- a.
  B[3,2] <- beta[1]
  B[4,2] <- beta[2]
  B[5,3] <- beta[3]
  B[5,4] <- beta[4]
  B[6,5] <- b
  
  I <- diag(6)
  evec <- diag(B %*% t(B))
  evec[5] <- evec[5] + 2*beta[1]*beta[2]*beta[3]*beta[4]
  
  Psi <- diag(6) - diag(evec)
  Sigma <- solve(I - B) %*% Psi %*% t(solve(I - B))
  #order = x, m1, m2, m3, m4, y
  
  #add a column for M (sum of m1-m4)
  sumMat <- matrix(0, 6, 7)
  sumMat[1:6,1:6] <- diag(6)
  sumMat[2:5,7] <- 1
  
  Sigma <- t(sumMat) %*% Sigma %*% sumMat
  #order = x, m1, m2, m3, m4, y, M
  
  #create unreliable indicators of m1-m4, call them y1-y4:
  Lambda <- matrix(0, 11, 7)
  diag(Lambda[1:7, 1:7]) <- 1
  diag(Lambda[8:11, 1:4]) <- sqrt(r)
  Theta <- matrix(0, 11, 11)
  diag(Theta[8:11, 8:11]) <- 1-r
  
  Sigma2 <- Lambda %*% Sigma %*% t(Lambda) + Theta
  #order = x, m1, m2, m3, m4, y, M, y1, y2, y3, y4
  
  #add a column for M (sum of y1-y4)
  sumMat <- matrix(0, 11, 12)
  sumMat[1:11,1:11] <- diag(11)
  sumMat[8:11,12] <- 1
  
  Sigma3 <- t(sumMat) %*% Sigma2 %*% sumMat  
  #order = x, m1, m2, m3, m4, y, M, y1, y2, y3, y4, M2
  
  colnames(Sigma3) <- c("X", "m1", "m2", "m3", "m4", "Y", "M", "y1", "y2", "y3", "y4", "M2")
  
  return(Sigma3)    
}
################################################################################################### 