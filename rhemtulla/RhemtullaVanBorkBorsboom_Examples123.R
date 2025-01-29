##################################################################################### 
## R code to produce Examples 1-3 in Rhemtulla, M., van Bork, R., & Borsboom, D.  
## (under review at Psych Methods)
## Written by Mijke Rhemtulla
## email: mrhemtulla@ucdavis.edu if you have comments/questions
#####################################################################################


# preliminary stuff #################################################################
source("RhemtullaVanBorkBorsboom_Functions.R")
windowsFonts(A=windowsFont("Verdana"))
#####################################################################################

## Example 1: causal indicators with measurement error ########################
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

# generate a single pop cov matrix and run models: 
a <- .6
b <- .6
res.M <- .25
avec <- c(.35, .45, .55, .65)
r <- .8
data <- generateCov3(a = a, b = b, avec = avec, covm = NA, res.M = res.M, r = .8)  

meana <- mean(avec)
vara <- var(avec)
rangea <- range(avec)[2] - range(avec)[1]
amat <- avec %*% t(avec)
gamma <- a/sum(avec)
covm <- (1 - 4*gamma^2 - res.M)/(12*gamma^2)

#get parameter values and model results
avec
covm
gamma

evec <- 1 - avec^2
emat <- covm - amat #residual covariances among m1-m4: 

fitLVMeasurementMod <- cfa (model = 'LM =~ m1 + m2 + m3 + m4', sample.cov = data, 
                            sample.nobs = 100000)

fitCompA <- sem(model = CompModA, sample.cov = data, sample.nobs = 100000)  
fitLVA <- sem(model = LVModA, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

fitCompB <- sem(model = CompModB, sample.cov = data, sample.nobs = 100000)  
fitLVB <- sem(model = LVModB, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

fitCompMed <- sem(model = CompModMed, sample.cov = data, sample.nobs = 100000)  
fitLVMed <- sem(model = LVModMed, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

#LV Measurement Model (CFA)
inspect(fitLVMeasurementMod, "std")$lambda[1:4,1]
sqrt((inspect(fitLVMeasurementMod, "fit")["fmin"]*2)/inspect(fitLVMeasurementMod,"fit")["df"])#pop RMSEA

#A-only models: 
inspect(fitLVA, "std")$lambda[1:4,1] #factor loadings
inspect(fitCompA, "std")$beta[1,2] #a path A-only composite model  
inspect(fitLVA, "std")$beta[1,2] #a path A-only LV model
sqrt((inspect(fitCompA, "fit")["fmin"]*2)/inspect(fitCompA,"fit")["df"])#pop RMSEA
sqrt((inspect(fitLVA, "fit")["fmin"]*2)/inspect(fitLVA,"fit")["df"])#pop RMSEA


#B-only models:
inspect(fitLVB, "std")$lambda[1:4,1] #factor loadings
inspect(fitCompB, "std")$beta[1,2] #b path B-only composite model 
inspect(fitLVB, "std")$beta[2,1] #b path B-only LV model 
sqrt((inspect(fitLVB, "fit")["fmin"]*2)/inspect(fitLVB,"fit")["df"])#pop RMSEA
sqrt((inspect(fitCompB, "fit")["fmin"]*2)/inspect(fitCompB,"fit")["df"])#pop RMSEA


#full mediation models: 
inspect(fitLVMed, "std")$lambda[1:4,1]
inspect(fitCompMed, "std")$beta[1,3] #a path med composite model 
inspect(fitLVMed, "std")$beta[1,3] #a path med LV model
inspect(fitCompMed, "std")$beta[2,1] #b path med composite model  
inspect(fitLVMed, "std")$beta[2,1] #b path med LV model 
parameterEstimates(fitCompMed, standardized = TRUE)[6,12] #ab path med composite model
parameterEstimates(fitLVMed, standardized = TRUE)[14,12] #ab path med LV model std
sqrt((inspect(fitLVMed, "fit")["fmin"]*2)/inspect(fitLVMed,"fit")["df"])#pop RMSEA
sqrt((inspect(fitCompMed, "fit")["fmin"]*2)/inspect(fitCompMed,"fit")["df"])#pop RMSEA
#####################################################################################

## Example 2: MIMIC model with measurement error ####################################
# generate a single pop cov matrix and run models: 

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

#why does every value return weird values of covm and negative residual cov? 
av <- seq(0, 1, by = .01) #average of avec
ga <- .6/(2*av)            #gamma
ga2 <- .4/(2*av)           #gamma if a = .4
plot(av, ga, type = "l", ylim = c(0, 1))
lines(av, ga2, col = "grey")
cm <- (.75 - 2*ga^2)/(2*ga^2)     #covm if a = .6
cm2 <- (.75 - 2*ga2^2)/(2*ga2^2)  #covm if a = .4
lines(av, cm, col = "red")
lines(av, cm2, col = "pink")
rc <- cm - av^2 #residual cov
rc2 <- cm2 - av^2 #residual cov

lines(av, rc, col = "darkgreen")
lines(av, rc2, col = "lightgreen")

#so to get positive values of everything, mean avec should be >= .6 or >= .4 if a = .4

a <- .6
b <- .6
res.M <- .25
r <- .8
avec <- c(.55, .65)
gamma <- a/sum(avec)
covm <- (1 - 2*gamma^2 - res.M)/(2*gamma^2) 
covm
evec <- 1 - avec^2
emat <- covm - amat 

lamvec <- c(.55, .65) 

data <- generateCov.mimic.r(a = a, b = b, avec = avec, lamvec = lamvec, res.M = res.M, covm = NA, r = .8) 
data[c(8:9,3:4),c(8:9,3:4)]

fitLVMeasurementMod <- cfa (model = 'LM =~ y1 + y2 + m3 + m4', sample.cov = data, 
                            sample.nobs = 100000)

fitCompA <- sem(model = CompModA, sample.cov = data, sample.nobs = 100000)  
fitLVA <- sem(model = LVModA, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

fitCompB <- sem(model = CompModB, sample.cov = data, sample.nobs = 100000)  
fitLVB <- sem(model = LVModB, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

fitCompMed <- sem(model = CompModMed, sample.cov = data, sample.nobs = 100000)  
fitLVMed <- sem(model = LVModMed, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

#LV Measurement Model (CFA)
inspect(fitLVMeasurementMod, "std")$lambda[1:4,1]
sqrt((inspect(fitLVMeasurementMod, "fit")["fmin"]*2)/inspect(fitLVMeasurementMod,"fit")["df"])#pop RMSEA

#A-only models: 
inspect(fitLVA, "std")$lambda[1:4,1] #factor loadings
inspect(fitCompA, "std")$beta[1,2] #a path A-only composite model  
inspect(fitLVA, "std")$beta[1,2] #a path A-only LV model
sqrt((inspect(fitCompA, "fit")["fmin"]*2)/inspect(fitCompA,"fit")["df"])#pop RMSEA
sqrt((inspect(fitLVA, "fit")["fmin"]*2)/inspect(fitLVA,"fit")["df"])#pop RMSEA

#B-only models:
inspect(fitLVB, "std")$lambda[1:4,1] #factor loadings
inspect(fitCompB, "std")$beta[1,2] #b path B-only composite model 
inspect(fitLVB, "std")$beta[2,1] #b path B-only LV model 
sqrt((inspect(fitLVB, "fit")["fmin"]*2)/inspect(fitLVB,"fit")["df"])#pop RMSEA
sqrt((inspect(fitCompB, "fit")["fmin"]*2)/inspect(fitCompB,"fit")["df"])#pop RMSEA

#full mediation models: 
inspect(fitLVMed, "std")$lambda[1:4,1]
inspect(fitCompMed, "std")$beta[1,3] #a path med composite model 
inspect(fitLVMed, "std")$beta[1,3] #a path med LV model
inspect(fitCompMed, "std")$beta[2,1] #b path med composite model  
inspect(fitLVMed, "std")$beta[2,1] #b path med LV model 
parameterEstimates(fitCompMed, standardized = TRUE)[6,12] #ab path med composite model
parameterEstimates(fitLVMed, standardized = TRUE)[14,12] #ab path med LV model std
sqrt((inspect(fitLVMed, "fit")["fmin"]*2)/inspect(fitLVMed,"fit")["df"])#pop RMSEA
sqrt((inspect(fitCompMed, "fit")["fmin"]*2)/inspect(fitCompMed,"fit")["df"])#pop RMSEA
#####################################################################################

## Example 3: network indicators with measurement error #############################

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


# generate a single pop cov matrix and run models: 
a <- .6  
b <- .6
beta <- c(.5, .6, .6, .7)
data <- generateCov.net.r(a = a, b = b, beta = beta, r = .8) 
eigen(data[-c(7:12),-c(7:12)])$values
a. <- data[1,2]
a. #.833

data[8:11,8:11]

fitLVMeasurementMod <- cfa (model = 'LM =~ y1 + y2 + y3 + y4', sample.cov = data, 
                            sample.nobs = 100000)

fitCompA <- sem(model = CompModA, sample.cov = data, sample.nobs = 100000)  
fitLVA <- sem(model = LVModA, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

fitCompB <- sem(model = CompModB, sample.cov = data, sample.nobs = 100000)  
fitLVB <- sem(model = LVModB, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

fitCompMed <- sem(model = CompModMed, sample.cov = data, sample.nobs = 100000)  
fitLVMed <- sem(model = LVModMed, sample.cov = data, sample.nobs = 100000, std.lv = TRUE)  

#LV Measurement Model (CFA)
inspect(fitLVMeasurementMod, "std")$lambda[1:4,1]
sqrt((inspect(fitLVMeasurementMod, "fit")["fmin"]*2)/inspect(fitLVMeasurementMod,"fit")["df"])#pop RMSEA

#A-only models: 
inspect(fitLVA, "std")$lambda[1:4,1] #factor loadings
inspect(fitCompA, "std")$beta[1,2] #a path A-only composite model  
inspect(fitLVA, "std")$beta[1,2] #a path A-only LV model
sqrt((inspect(fitCompA, "fit")["fmin"]*2)/inspect(fitCompA,"fit")["df"])#pop RMSEA
sqrt((inspect(fitLVA, "fit")["fmin"]*2)/inspect(fitLVA,"fit")["df"])#pop RMSEA

#B-only models:
inspect(fitLVB, "std")$lambda[1:4,1] #factor loadings
inspect(fitCompB, "std")$beta[1,2] #b path B-only composite model 
inspect(fitLVB, "std")$beta[2,1] #b path B-only LV model 
sqrt((inspect(fitLVB, "fit")["fmin"]*2)/inspect(fitLVB,"fit")["df"])#pop RMSEA
sqrt((inspect(fitCompB, "fit")["fmin"]*2)/inspect(fitCompB,"fit")["df"])#pop RMSEA

#full mediation models: 
inspect(fitLVMed, "std")$lambda[1:4,1]
inspect(fitCompMed, "std")$beta[1,3] #a path med composite model 
inspect(fitLVMed, "std")$beta[1,3] #a path med LV model
inspect(fitCompMed, "std")$beta[2,1] #b path med composite model  
inspect(fitLVMed, "std")$beta[2,1] #b path med LV model 
parameterEstimates(fitCompMed, standardized = TRUE)[6,12] #ab path med composite model
parameterEstimates(fitLVMed, standardized = TRUE)[14,12] #ab path med LV model std
sqrt((inspect(fitLVMed, "fit")["fmin"]*2)/inspect(fitLVMed,"fit")["df"])#pop RMSEA
sqrt((inspect(fitCompMed, "fit")["fmin"]*2)/inspect(fitCompMed,"fit")["df"])#pop RMSEA
#####################################################################################