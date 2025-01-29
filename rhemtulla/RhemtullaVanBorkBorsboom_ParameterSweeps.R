############################################################################### 
## R code to produce Figures for Rhemtulla, M., van Bork, R., & Borsboom, D.  
## (under review at Psych Methods)
## Written by Mijke Rhemtulla
## email: mrhemtulla@ucdavis.edu if you have comments/questions
############################################################################### 

#Example 1: causal indicator model with unreliable indicators and residual variance on M

# preliminary stuff #####################################################################################
source("RhemtullaVanBorkBorsboom_Functions.R")
windowsFonts(A=windowsFont("Verdana"))
######################################################################################################

# Example 1 Parameter Sweep ##########################################################################

#models to fit: 
CompModMed <- 'M2 ~ a*X
Y ~ b*M2
ab := a*b'

LVModMed <- 'LM =~ y1 + y2 + y3 + y4
LM ~ a*X
Y ~ b*LM
ab := a*b'

CompModA <- 'M2 ~ A*X'
LVModA <- 'LM =~ y1 + y2 + y3 + y4
LM ~ a*X'

CompModB <- 'Y ~ b*M2'
LVModB <- 'LM =~ y1 + y2 + y3 + y4
Y ~ b*LM'


# fx to sample random values of a then generate population covariance matrix, fit models, save results: 
simEx1 <- function(i){
  a <- .6
  b <- .6
  res.M <- .25
  avec <- runif(4, 0, .8) #in sim2 this had max = .5
  data <- generateCov3(a = a, b = b, avec = avec, covm = NA, res.M = res.M, r = .8)  
  while (min(eigen(data[-c(5,8),-c(5,8)])$values) < 1e-10 || data[1,2] < 0) { #ensure that sigma is positive definite! 
    avec <- runif(4, 0, .8) #in sim2 this had max = .5
    data <- generateCov3(a = a, b = b, avec = avec, covm = NA, res.M = res.M, r = .8)
  }
  SEs.LVA <- TRUE  
  SEs.LVB <- TRUE  
  SEs.LVMed <- TRUE  

  fitCompA <- withTimeout(sem(model = CompModA, sample.cov = data, sample.nobs = 100000),
                          timeout = 10, onTimeout = "silent")
  fitLVA <- withTimeout(sem(model = LVModA, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                        timeout = 10, onTimeout = "silent")
  conv.LVA <- try(inspect(fitLVA, "converged"))
  if (conv.LVA == FALSE) {return(rep(NA, 37))}
  if (class(try(vcov(fitLVA))) == "try-error") {return(rep(NA, 37))}
  
  fitCompB <- withTimeout(sem(model = CompModB, sample.cov = data, sample.nobs = 100000),
                          timeout = 10, onTimeout = "silent")
  fitLVB <- withTimeout(sem(model = LVModB, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                        timeout = 10, onTimeout = "silent")
  conv.LVB <- try(inspect(fitLVB, "converged"))
  if (conv.LVB == FALSE) {return(rep(NA, 37))}
  if (class(try(vcov(fitLVB))) == "try-error") {return(rep(NA, 37))}
  
  fitCompMed <- withTimeout(sem(model = CompModMed, sample.cov = data, sample.nobs = 100000),
                            timeout = 10, onTimeout = "silent")
  fitLVMed <- withTimeout(sem(model = LVModMed, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                          timeout = 10, onTimeout = "silent")
  conv.LVMed <- try(inspect(fitLVMed, "converged"))
  if (conv.LVMed == FALSE) {return(rep(NA, 37))}
  if (class(try(vcov(fitLVMed))) == "try-error") {return(rep(NA, 37))}
  
  fit.LVA <- withTimeout(inspect(fitLVA, "fit")["fmin"],
                         timeout = 10, onTimeout = "silent")
  fit.LVB <- withTimeout(inspect(fitLVB, "fit")["fmin"],
                         timeout = 10, onTimeout = "silent")
  fit.LVMed <- withTimeout(inspect(fitLVMed, "fit")["fmin"],
                           timeout = 10, onTimeout = "silent")
  
  meana <- mean(avec)
  vara <- var(avec)
  rangea <- range(avec)[2] - range(avec)[1]
  amat <- avec %*% t(avec)
  gamma <- a/sum(avec)
  covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)   #NEW so total var(M) includes residual variance
  
  return(c(meana, rangea, vara, covm, gamma,
           
           inspect(fitLVA, "std")$lambda[1:4,1],
           inspect(fitCompA, "std")$beta[1,2], #a path A-only composite model  
           inspect(fitLVA, "std")$beta[1,2], #a path A-only LV model
           
           inspect(fitLVB, "std")$lambda[1:4,1],
           inspect(fitCompB, "std")$beta[1,2], #b path B-only composite model 
           inspect(fitLVB, "std")$beta[2,1], #b path B-only LV model 
           
           inspect(fitCompMed, "std")$beta[1,3], #a path med composite model 
           inspect(fitCompMed, "std")$beta[2,1], #b path med composite model  
           parameterEstimates(fitCompMed, standardized = TRUE)[6,12], #ab path med composite model
           
           inspect(fitLVMed, "std")$lambda[1:4,1],  
           inspect(fitLVMed, "std")$beta[1,3], #a path med LV model
           inspect(fitLVMed, "std")$beta[2,1], #b path med LV model 
           parameterEstimates(fitLVMed, standardized = TRUE)[14,12], #ab path med LV model std
           
           inspect(fitCompMed, "fit")["fmin"], #fit of composite model
           fit.LVA, fit.LVB, fit.LVMed,
           conv.LVA, conv.LVB, conv.LVMed, SEs.LVA, SEs.LVB, SEs.LVMed)) 
}

#do this a lot of times: 
resultsEx1 <- matrix(NA, 50000, 37)
for (i in 1:40000){
  resultsEx1[i,] <- simEx1(i)
  print(i)
  flush.console()
}
   
results <- as.data.frame(resultsEx1)

#save(results, file = "Example1Sweep.Rda")
#load(file = "Example1Sweep.Rda")

colnames(results) <- c("mean.a", "range.a", "var.a", "covm", "gamma", 
                       "l1A", "l2A", "l3A", "l4A", "acompA", "aLVA",
                       "l1B", "l2B", "l3B", "l4B", "bcompB", "bLVB",
                       "acompM", "bcompM", "abcompM",
                       "l1M", "l2M", "l3M", "l4M", "aLVM", "bLVM", "abLVM", 
                       "fcompAB", "fit.LVA", "fit.LVB", "fit.LVM",
                       "conv.LVA", "conv.LVB", "conv.LVM", 
                       "SEs.LVA", "SEs.LVB", "SEs.LVM")

#discard rows that did not converge: 
results2 <- results[complete.cases(results),]
dim(results2)

results2 <- results2[1:20000,] #keep just the first 20000 assuming we got more than that

#Create a function to generate a continuous color palette
redPal <- colorRampPalette(c('#920000FF','#FFBDBDFF'))
purPal <- colorRampPalette(c('#490092FF','#D9B3FFFF'))

results2$col.var.a <- redPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col2.var.a <- purPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]

#make plots
jpeg("Ex1Sweep_abpaths.jpg", width = 8, height = 8, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, .5, .5, 1))
plot(results2$covm, results2$acompA, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(a) ~ path), 
     ylim = c(0, 1), pch = ".", cex = 2, col = "#490092FF", 
     main = "predictor-only model")
points(results2$covm, results2$aLVA, col = results2$col.var.a, pch = ".", cex = 2)
abline(h = .6)

legend(x = .2, y = .4, c("composite model", "common factor model"), 
       fill = c("#490092FF", "darkred"), bty = "n")

par(fig = c(.5, 1, .5, 1), new = TRUE)
plot(results2$covm, results2$bcompB, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(b) ~ path), 
     ylim = c(0, 1), pch = ".", cex = 2, col = "#490092FF",
     main = "outcome-only model")
points(results2$covm, results2$bLVB, col = results2$col.var.a, pch = ".", cex = 2)
abline(h = .6)

par(fig = c(0, .5, 0, .5), new = TRUE)
plot(results2$covm, results2$acompM, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(a) ~ path), 
     ylim = c(0, 1), pch = ".", cex = 2, col = "#490092FF",
     main = "mediation model")
points(results2$covm, results2$aLVM, col = results2$col.var.a, pch = ".", cex = 2)
abline(h = .6)

par(fig = c(.5, 1, 0, .5), new = TRUE)
plot(results2$covm, results2$bcompM, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(b) ~ path), 
     ylim = c(0, 1), pch = ".", cex = 2, col = "#490092FF",
     main = "mediation model")
points(results2$covm, results2$bLVM, col = results2$col.var.a, pch = ".", cex = 2)
abline(h = .6)

dev.off()

## plot model fit (RMSEA)
results2$RMSEA.LVM <- sqrt(results2$fit.LVM*2/9) 
results2$RMSEA.LVA <- sqrt(results2$fit.LVA*2/5) 
results2$RMSEA.LVB <- sqrt(results2$fit.LVB*2/5) 
results2$RMSEA.CompM <- sqrt(results2$fcompAB*2) 

jpeg("Ex1SweepRMSEA.jpg", width = 8, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, .5, 0, 1))
  plot(results2$covm, results2$RMSEA.LVA, type = "p", 
       xlab = "correlation between causal indicators",
       ylab = "Population RMSEA", ylim = c(0, 1), col = results2$col.var.a,  
       pch = ".", cex = 2, main = "predictor-only model")
  
  par(fig = c(.5, 1, 0, 1), new = TRUE)
  plot(results2$covm, results2$RMSEA.LVM, type = "p", 
       xlab = "correlation between causal indicators",
       ylab = "Population RMSEA", ylim = c(0, 1), col = results2$col.var.a, 
       pch = ".", cex = 2, main = "mediation model")
  points(results2$covm, results2$RMSEA.CompM, col = "#490092FF", pch = ".", cex = 2)
dev.off()

#plot reliability of unreliable observed composite as covm increases: 
r <- .8
covm <- seq(.2, .9, by = .001)
relObs <- (4*r + 12*r*covm)/(4 + 12*.8*covm) 
plot(covm, relObs, type = "l")

#correlation of a and b estimates in the lv med model: 
cor(results2$aLVM, results2$bLVM)

#correlation of mean factor loading and bias in the a path of pred-only model: 
results2$meanLam <- apply(results2[,6:9], 1, mean)
cor(results2$meanLam, results2$aLVA)

#correlation of bias in the a path of pred-only model and rmsea: 
cor(results2$RMSEA.LVA, results2$aLVA)
plot(results2$RMSEA.LVA, results2$aLVA, pch = ".")

######################################################################################################

# Example 2 Parameter Sweep ##########################################################################
CompModMed <- 'M2 ~ a*X
Y ~ b*M2
ab := a*b'

LVModMed <- 'LM =~ y1 + y2 + m3 + m4
LM ~ a*X
Y ~ b*LM
ab := a*b'

CompModA <- 'M2 ~ A*X'
LVModA <- 'LM =~ y1 + y2 + m3 + m4
LM ~ a*X'

CompModB <- 'Y ~ b*M2'
LVModB <- 'LM =~ y1 + y2 + m3 + m4
Y ~ b*LM'

simEx2 <- function(i){
  a <- .6
  b <- .6
  res.M <- .25
  r <- .8
  data <- matrix(1, 7, 7)
  while (min(eigen(data[-c(5,8),-c(5,8)])$values) < 1e-10 || data[1,2] < 0) {
    avec <- runif(2, .3, .8) 
    lamvec <- runif(2, .3, 1) 
    data <- generateCov.mimic.r(a = a, b = b, avec = avec, lamvec = lamvec, res.M = res.M, covm = NA, r = .8) 
  }
  SEs.LVA <- TRUE  
  SEs.LVB <- TRUE  
  SEs.LVMed <- TRUE  
  
  fitCompA <- withTimeout(sem(model = CompModA, sample.cov = data, sample.nobs = 100000),
                          timeout = 10, onTimeout = "silent")
  fitLVA <- withTimeout(sem(model = LVModA, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                        timeout = 10, onTimeout = "silent")
  conv.LVA <- try(inspect(fitLVA, "converged"))
  if (conv.LVA == FALSE) {return(rep(NA, 38))}
  if (class(try(vcov(fitLVA))) == "try-error") {return(rep(NA, 38))}
  
  fitCompB <- withTimeout(sem(model = CompModB, sample.cov = data, sample.nobs = 100000),
                          timeout = 10, onTimeout = "silent")
  fitLVB <- withTimeout(sem(model = LVModB, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                        timeout = 10, onTimeout = "silent")
  conv.LVB <- try(inspect(fitLVB, "converged"))
  if (conv.LVB == FALSE) {return(rep(NA, 38))}
  if (class(try(vcov(fitLVB))) == "try-error") {return(rep(NA, 38))}
  
  fitCompMed <- withTimeout(sem(model = CompModMed, sample.cov = data, sample.nobs = 100000),
                            timeout = 10, onTimeout = "silent")
  conv.CompMed <- try(inspect(fitCompMed, "converged"))
  if (conv.CompMed == FALSE) {return(rep(NA, 38))}
  fitLVMed <- withTimeout(sem(model = LVModMed, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                          timeout = 10, onTimeout = "silent")
  conv.LVMed <- try(inspect(fitLVMed, "converged"))
  if (conv.LVMed == FALSE) {return(rep(NA, 38))}
  if (class(try(vcov(fitLVMed))) == "try-error") {return(rep(NA, 38))}
  
  fit.LVA <- withTimeout(inspect(fitLVA, "fit")["fmin"],
                         timeout = 10, onTimeout = "silent")
  fit.LVB <- withTimeout(inspect(fitLVB, "fit")["fmin"],
                         timeout = 10, onTimeout = "silent")
  fit.LVMed <- withTimeout(inspect(fitLVMed, "fit")["fmin"],
                           timeout = 10, onTimeout = "silent")
  gamma <- a/sum(avec)
  covm <- (1 - 2*gamma^2 - res.M)/(2*gamma^2)   
  
  return(c(avec, lamvec, gamma, covm,
           
           inspect(fitLVA, "std")$lambda[1:4,1],
           inspect(fitCompA, "std")$beta[1,2], #a path A-only composite model  
           inspect(fitLVA, "std")$beta[1,2], #a path A-only LV model
           
           inspect(fitLVB, "std")$lambda[1:4,1],
           inspect(fitCompB, "std")$beta[1,2], #b path B-only composite model 
           inspect(fitLVB, "std")$beta[2,1], #b path B-only LV model 
           
           inspect(fitCompMed, "std")$beta[1,3], #a path med composite model 
           inspect(fitCompMed, "std")$beta[2,1], #b path med composite model  
           parameterEstimates(fitCompMed, standardized = TRUE)[6,12], #ab path med composite model
           
           inspect(fitLVMed, "std")$lambda[1:4,1],  
           inspect(fitLVMed, "std")$beta[1,3], #a path med LV model
           inspect(fitLVMed, "std")$beta[2,1], #b path med LV model 
           parameterEstimates(fitLVMed, standardized = TRUE)[14,12], #ab path med LV model std
           
           inspect(fitCompMed, "fit")["fmin"], #fit of composite model
           fit.LVA, fit.LVB, fit.LVMed,
           conv.LVA, conv.LVB, conv.LVMed, SEs.LVA, SEs.LVB, SEs.LVMed)) 
}

resultsEx2 <- matrix(NA, 25000, 38)
for (i in 1:25000){
  resultsEx2[i,] <- simEx2(i)
  print(i)
  flush.console()
}

results <- as.data.frame(resultsEx2)

#save(results, file = "Example2Sweep.Rda")
#load(file = "Example2Sweep.Rda")

colnames(results) <- c("a1", "a2", "lam1", "lam2", "gamma", "covm", 
                         "l1A", "l2A", "l3A", "l4A", "acompA", "aLVA",
                         "l1B", "l2B", "l3B", "l4B", "bcompB", "bLVB",
                         "acompM", "bcompM", "abcompM",
                         "l1M", "l2M", "l3M", "l4M", "aLVM", "bLVM", "abLVM", 
                         "fcompAB", "fit.LVA", "fit.LVB", "fit.LVM",
                         "conv.LVA", "conv.LVB", "conv.LVM", 
                         "SEs.LVA", "SEs.LVB", "SEs.LVM")


#discard rows that did not converge: 
results2 <- results[complete.cases(results),]
dim(results2)
results2 <- results2[1:20000,]

#Create a function to generate a continuous color palette
redPal <- colorRampPalette(c('#920000FF','#FFBDBDFF'))
purPal <- colorRampPalette(c('#490092FF','#D9B3FFFF'))

results2$var.a <- apply(cbind(results2$a1, results2$a2), 1, var)
results2$col.var.a <- redPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]
results2$col2.var.a <- purPal(10)[as.numeric(cut(results2$var.a,breaks = 10))]

#make plots

jpeg("Ex2Sweep_abpaths.jpg", width = 8, height = 8, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, .5, .5, 1))
plot(results2$covm, results2$acompA, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(a) ~ path), 
     col = results2$col2.var.a, ylim = c(0, 1), pch = ".", cex = 1,
     main = "predictor-only model", xlim = c(0, 1))
points(results2$covm, results2$aLVA, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .6)

legend(x = 0, y = .4, c("composite model", "common factor model"), 
       fill = c("#490092FF", "darkred"), bty = "n")

par(fig = c(.5, 1, .5, 1), new = TRUE)
plot(results2$covm, results2$bcompB, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(b) ~ path), 
     col = results2$col2.var.a, ylim = c(0, 1), pch = ".", cex = 1,
     main = "outcome-only model", xlim = c(0, 1))
points(results2$covm, results2$bLVB, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .6)

par(fig = c(0, .5, 0, .5), new = TRUE)
plot(results2$covm, results2$acompM, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(a) ~ path), 
     col = results2$col2.var.a, ylim = c(0, 1), pch = ".", cex = 1,
     main = "mediation model", xlim = c(0, 1))
points(results2$covm, results2$aLVM, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .6)

par(fig = c(.5, 1, 0, .5), new = TRUE)
plot(results2$covm, results2$bcompM, type = "p", 
     xlab = expression(correlation ~ between ~ causal ~ indicators),
     ylab = expression(estimated ~ italic(b) ~ path), 
     col = results2$col2.var.a, ylim = c(0, 1), pch = ".", cex = 1,
     main = "mediation model", xlim = c(0, 1))
points(results2$covm, results2$bLVM, col = results2$col.var.a, pch = ".", cex = 1)
abline(h = .6)
dev.off()

## plot model fit (RMSEA)
results2$RMSEA.LVM <- sqrt(results2$fit.LVM*2/9) 
results2$RMSEA.LVA <- sqrt(results2$fit.LVA*2/5) 
results2$RMSEA.LVB <- sqrt(results2$fit.LVB*2/5) 
results2$RMSEA.CompM <- sqrt(results2$fcompAB*2) 

jpeg("Ex2SweepRMSEA.jpg", width = 12, height = 4, units = "in", res = 800)
  par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
  par(fig = c(0, 1/3, 0, 1))
  plot(results2$covm, results2$RMSEA.LVA, type = "p", 
       xlab = "correlation between causal indicators",
       ylab = "Population RMSEA", col = results2$col.var.a, ylim = c(0, 1), 
       pch = ".", cex = 1, main = "predictor-only model", xlim = c(0, 1))
  points(results2$covm, results2$RMSEA.CompA, col = results2$col2.var.a, pch = ".")
  
  par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
  plot(results2$covm, results2$RMSEA.LVB, type = "p", 
       xlab = "correlation between causal indicators",
       ylab = "Population RMSEA", col = results2$col.var.a, ylim = c(0, 1), 
       pch = ".", cex = 1, main = "outcome-only model", xlim = c(0, 1))
  points(results2$covm, results2$RMSEA.CompB, col = results2$col2.var.a, pch = ".")
  
  par(fig = c(2/3, 1, 0, 1), new = TRUE)
  plot(results2$covm, results2$RMSEA.LVM, type = "p", 
       xlab = "correlation between causal indicators",
       ylab = "Population RMSEA", col = results2$col.var.a, ylim = c(0, 1), 
       pch = ".", cex = 1, main = "mediation model", xlim = c(0, 1))
  points(results2$covm, results2$RMSEA.CompM, col = results2$col2.var.a, pch = ".")

dev.off()

#correlation of a and b estimates in the lv med model: 
cor(results2$aLVM, results2$bLVM)

#correlation of mean factor loading and bias in the a path of pred-only model: 
results2$meanLam <- apply(results2[,6:9], 1, mean)
cor(results2$meanLam, results2$aLVA)

#correlation of bias in the a path of pred-only model and rmsea: 
cor(results2$RMSEA.LVA, results2$aLVA)
plot(results2$RMSEA.LVA, results2$aLVA, pch = ".")

######################################################################################################

# Example 3 Parameter Sweep ##########################################################################
simEx3 <- function(i){
  a <- .6  
  b <- .6
  
  data <- matrix(1, 7, 7)
  while(min(eigen(data[-7,-7])$values) <= 0){
    beta <- runif(4, 0, .8) 
    data <- generateCov.net.r(a = a, b = b, beta = beta, r = .8) 
  }
  
  a. <- data[1,2]
  
  SEs.LVA <- TRUE  
  SEs.LVB <- TRUE  
  SEs.LVMed <- TRUE  
  
  fitCompA <- withTimeout(sem(model = CompModA, sample.cov = data, sample.nobs = 100000),
                          timeout = 10, onTimeout = "silent")
  fitLVA <- withTimeout(sem(model = LVModA, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                        timeout = 10, onTimeout = "silent")
  conv.LVA <- try(inspect(fitLVA, "converged"))
  if (conv.LVA == FALSE) {return(rep(NA, 37))}
  if (class(try(vcov(fitLVA))) == "try-error") {return(rep(NA, 37))}
  
  fitCompB <- withTimeout(sem(model = CompModB, sample.cov = data, sample.nobs = 100000),
                          timeout = 10, onTimeout = "silent")
  fitLVB <- withTimeout(sem(model = LVModB, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                        timeout = 10, onTimeout = "silent")
  conv.LVB <- try(inspect(fitLVB, "converged"))
  if (conv.LVB == FALSE) {return(rep(NA, 37))}
  if (class(try(vcov(fitLVB))) == "try-error") {return(rep(NA, 37))}
  
  fitCompMed <- withTimeout(sem(model = CompModMed, sample.cov = data, sample.nobs = 100000),
                            timeout = 10, onTimeout = "silent")
  conv.CompMed <- try(inspect(fitCompMed, "converged"))
  if (conv.CompMed == FALSE) {return(rep(NA, 37))}
  fitLVMed <- withTimeout(sem(model = LVModMed, sample.cov = data, sample.nobs = 100000, std.lv = TRUE),
                          timeout = 10, onTimeout = "silent")
  conv.LVMed <- try(inspect(fitLVMed, "converged"))
  if (conv.LVMed == FALSE) {return(rep(NA, 37))}
  if (class(try(vcov(fitLVMed))) == "try-error") {return(rep(NA, 37))}
  
  fit.LVA <- withTimeout(inspect(fitLVA, "fit")["fmin"],
                         timeout = 10, onTimeout = "silent")
  fit.LVB <- withTimeout(inspect(fitLVB, "fit")["fmin"],
                         timeout = 10, onTimeout = "silent")
  fit.LVMed <- withTimeout(inspect(fitLVMed, "fit")["fmin"],
                           timeout = 10, onTimeout = "silent")
  return(c(a., beta,  
           
           inspect(fitLVA, "std")$lambda[1:4,1],
           inspect(fitCompA, "std")$beta[1,2], #a path A-only composite model  
           inspect(fitLVA, "std")$beta[1,2], #a path A-only LV model
           
           inspect(fitLVB, "std")$lambda[1:4,1],
           inspect(fitCompB, "std")$beta[1,2], #b path B-only composite model 
           inspect(fitLVB, "std")$beta[2,1], #b path B-only LV model 
           
           inspect(fitCompMed, "std")$beta[1,3], #a path med composite model 
           inspect(fitCompMed, "std")$beta[2,1], #b path med composite model  
           parameterEstimates(fitCompMed, standardized = TRUE)[6,12], #ab path med composite model
           
           inspect(fitLVMed, "std")$lambda[1:4,1],  
           inspect(fitLVMed, "std")$beta[1,3], #a path med LV model
           inspect(fitLVMed, "std")$beta[2,1], #b path med LV model 
           parameterEstimates(fitLVMed, standardized = TRUE)[14,12], #ab path med LV model std
           
           inspect(fitCompMed, "fit")["fmin"], #fit of composite model
           fit.LVA, fit.LVB, fit.LVMed,
           conv.LVA, conv.LVB, conv.LVMed, SEs.LVA, SEs.LVB, SEs.LVMed)) 
}

CompModMed <- 'M2 ~ a*X
Y ~ b*M2
ab := a*b'

LVModMed <- 'LM =~ y1 + y2 + y3 + y4
LM ~ a*X
Y ~ b*LM
ab := a*b'

CompModA <- 'M2 ~ A*X'
LVModA <- 'LM =~ y1 + y2 + y3 + y4
LM ~ a*X'

CompModB <- 'Y ~ b*M2'
LVModB <- 'LM =~ y1 + y2 + y3 + y4
Y ~ b*LM'

#resultsEx3 <- matrix(NA, 50000, 37)
for (i in 450:25000){
  resultsEx3[i,] <- simEx3(i)
  print(i)
  flush.console()
}

#results <- as.data.frame(resultsEx3)
#save(results, file = "Example3Sweep.Rda")
#load(file = "Example3Sweep.Rda")

colnames(results) <- c("a.", "b1", "b2", "b3", "b4", 
                         "l1A", "l2A", "l3A", "l4A", "acompA", "aLVA",
                         "l1B", "l2B", "l3B", "l4B", "bcompB", "bLVB",
                         "acompM", "bcompM", "abcompM",
                         "l1M", "l2M", "l3M", "l4M", "aLVM", "bLVM", "abLVM", 
                         "fcompAB", "fit.LVA", "fit.LVB", "fit.LVM",
                         "conv.LVA", "conv.LVB", "conv.LVM", 
                         "SEs.LVA", "SEs.LVB", "SEs.LVM")
results3a <- data.frame(results)


#discard rows that did not converge: 
results2 <- results[complete.cases(results),]
dim(results2) #24922
results2 <- results2[1:20000,]

#Create a function to generate a continuous color palette
redPal <- colorRampPalette(c('#920000FF','#FFBDBDFF'))
purPal <- colorRampPalette(c('#490092FF','#D9B3FFFF'))

#make plots

jpeg("Ex3Sweep_abpaths.jpg", width = 12, height = 8, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, .5, 1))
plot(results2$a., results2$acompA, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = expression(estimated ~ italic(a) ~ path), 
     col = "#490092FF", ylim = c(0, 1), pch = ".", cex = 1,
     main = "predictor-only model", xlim = c(.7, 1))
points(results2$a., results2$aLVA, col = "#920000FF", pch = ".", cex = 1)
#abline(h = .6)

legend(x = .70, y = .5, c("composite model", "common factor model"), 
       fill = c("#490092FF", "#920000FF"), bty = "n")

par(fig = c(1/3, 2/3, .5, 1), new = TRUE)
plot(results2$a., results2$bcompB, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = expression(estimated ~ italic(b) ~ path), 
     col = "#490092FF", ylim = c(0, 1), pch = ".", cex = 1,
     main = "outcome-only model", xlim = c(.7, 1))
points(results2$a., results2$bLVB, col = "#920000FF", pch = ".", cex = 1)
#abline(h = .6)

par(fig = c(0, 1/3, 0, .5), new = TRUE)
plot(results2$a., results2$acompM, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = expression(estimated ~ italic(a) ~ path), 
     col = "#490092FF", ylim = c(0, 1), pch = ".", cex = 1,
     main = "mediation model", xlim = c(.7, 1))
points(results2$a., results2$aLVM, col = "#920000FF", pch = ".", cex = 1)
#abline(h = .6)

par(fig = c(1/3, 2/3, 0, .5), new = TRUE)
plot(results2$a., results2$bcompM, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = expression(estimated ~ italic(b) ~ path), 
     col = "#490092FF", ylim = c(0, 1), pch = ".", cex = 1,
     main = "mediation model", xlim = c(.7, 1))
points(results2$a., results2$bLVM, col = "#920000FF", pch = ".", cex = 1)
#abline(h = .6)

par(fig = c(2/3, 1, 0, .5), new = TRUE)
plot(results2$a., results2$abcompM, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = expression(estimated ~ italic(ab) ~ path), 
     col = "#490092FF", ylim = c(0, 1), pch = ".", cex = 1,
     main = "mediation model", xlim = c(.7, 1))
points(results2$a., results2$abLVM, col = "#920000FF", pch = ".", cex = 1)
abline (h = .36)
dev.off()

## plot model fit (RMSEA)
results2$RMSEA.LVM <- sqrt(results2$fit.LVM*2/9) 
results2$RMSEA.LVA <- sqrt(results2$fit.LVA*2/5) 
results2$RMSEA.LVB <- sqrt(results2$fit.LVB*2/5) 
results2$RMSEA.CompM <- sqrt(results2$fcompAB*2) 

jpeg("Ex3SweepRMSEA.jpg", width = 12, height = 4, units = "in", res = 800)
par(mar = c(5, 5, 1, 1), family = "A", mgp = c(3, 1, 0), xaxs = "i", yaxs = "i")
par(fig = c(0, 1/3, 0, 1))
plot(results2$a., results2$RMSEA.LVA, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = "Population RMSEA", 
     col = "#920000FF", ylim = c(0, 1), 
     pch = ".", cex = 1, main = "predictor-only model", xlim = c(.7, 1))
points(results2$a., results2$RMSEA.CompA, col = "#490092FF", pch = ".")

legend(x = .70, y = .85, c("composite model", "common factor model"), 
       fill = c("#490092FF", "#920000FF"), bty = "n")


par(fig = c(1/3, 2/3, 0, 1), new = TRUE)
plot(results2$a., results2$RMSEA.LVB, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = "Population RMSEA", col = "#920000FF", ylim = c(0, 1), 
     pch = ".", cex = 1, main = "outcome-only model", xlim = c(.7, 1))
points(results2$a1a4, results2$RMSEA.CompB, col = "#490092FF", pch = ".")

par(fig = c(2/3, 1, 0, 1), new = TRUE)
plot(results2$a., results2$RMSEA.LVM, type = "p", 
     xlab = expression(italic(X)%->%italic(m[1]) ~ path ~ value),
     ylab = "Population RMSEA", col = "#920000FF", ylim = c(0, 1), 
     pch = ".", cex = 1, main = "mediation model", xlim = c(.7, 1))
points(results2$a., results2$RMSEA.CompM, col = "#490092FF", pch = ".")

dev.off()


#correlation of a and b estimates in the lv med model: 
cor(results2$aLVM, results2$bLVM)

#correlation of mean factor loading and bias in the a path of pred-only model: 
results2$meanLam <- apply(results2[,6:9], 1, mean)
cor(results2$meanLam, results2$aLVA)

#correlation of bias in the a path of pred-only model and rmsea: 
cor(results2$RMSEA.LVA, results2$aLVA)
plot(results2$RMSEA.LVA, results2$aLVA, pch = ".")

######################################################################################################