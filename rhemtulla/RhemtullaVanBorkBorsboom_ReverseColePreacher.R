############################################################################### 
## R code to produce replicate & reverse effects of Cole & Preacher (2014) 
## for Rhemtulla, M., van Bork, R., & Borsboom, D.  
## (under review at Psych Methods)
## Written by Mijke Rhemtulla
## email: mrhemtulla@ucdavis.edu if you have comments/questions
############################################################################### 

# preliminary stuff ###########################################################
library(lavaan)
###############################################################################

#reproduce Cole & Preacher's Model in Figure 3C: ###############################
TrueMod <- 'J ~ .3*F + .3*G + .3*H
            F ~ .3*D
            G ~ .3*D + .3*E
            H ~ .3*E
            D ~~ .3*E
            J ~~ 0.68302*J
            F ~~ .91*F
            G ~~ .766*G
            H ~~ .91*H
            D ~~ D'

data1 <- simulateData(TrueMod, sample.nobs = 100000, empirical = TRUE)
cov1 <- round(cov(data1), 3)

FitMod <- 'J ~ f*F + g*G + h*H
           F ~ b*D
           G ~ c*D + d*E
           H ~ e*E
           D ~~ a*E'

fitTrue <- sem(FitMod, sample.cov = cov1, sample.nobs = 100000)
summary(fitTrue, std = TRUE)
################################################################################

# split each variable into 3 components  #######################################
covJ <- matrix(NA, 9, 9)
colnames(covJ) <- c(colnames(cov1), paste("J", 1:3, sep = ""))
covJ[1:6, 1:6] <- cov1

covm <- .5 #correlation among components of J
covJi <- matrix(covm, 3, 3)
diag(covJi) <- 1
covJi <- covJi / sum(covJi) #rescale so unit-weighted sum has variance 1
covJ[7:9, 7:9] <- covJi

covJ[1:6, 7:9] <- covJ[1:6, 1]*(1/3)
covJ[7:9, 1:6] <- matrix(covJ[1, 1:6]*(1/3), 3, 6, byrow = TRUE)
#WOOOO this appears to work. Do the same for other letters: 

covF <- covJ
colnames(covF) <- c(colnames(cov1), paste("F", 1:3, sep = ""))
covF[1:6, 7:9] <- covF[1:6, 2]*(1/3)
covF[7:9, 1:6] <- matrix(covF[2, 1:6]*(1/3), 3, 6, byrow = TRUE)

covG <- covJ
colnames(covG) <- c(colnames(cov1), paste("G", 1:3, sep = ""))
covG[1:6, 7:9] <- covG[1:6, 3]*(1/3)
covG[7:9, 1:6] <- matrix(covG[3, 1:6]*(1/3), 3, 6, byrow = TRUE)

covH <- covJ
colnames(covH) <- c(colnames(cov1), paste("H", 1:3, sep = ""))
covH[1:6, 7:9] <- covH[1:6, 4]*(1/3)
covH[7:9, 1:6] <- matrix(covH[4, 1:6]*(1/3), 3, 6, byrow = TRUE)

covD <- covJ
colnames(covD) <- c(colnames(cov1), paste("D", 1:3, sep = ""))
covD[1:6, 7:9] <- covD[1:6, 5]*(1/3)
covD[7:9, 1:6] <- matrix(covD[5, 1:6]*(1/3), 3, 6, byrow = TRUE)

covE <- covJ
colnames(covE) <- c(colnames(cov1), paste("E", 1:3, sep = ""))
covE[1:6, 7:9] <- covE[1:6, 6]*(1/3)
covE[7:9, 1:6] <- matrix(covE[6, 1:6]*(1/3), 3, 6, byrow = TRUE)
################################################################################

# Fit true model and models w/ each variable replaced with reflective mod ######
FitMod <- 'J ~ f*F + g*G + h*H
           F ~ b*D
           G ~ c*D + d*E
           H ~ e*E
           D ~~ a*E'

fitTrue <- sem(FitMod, sample.cov = cov1, sample.nobs = 100000)
summary(fitTrue, std = TRUE)
# everything accurate

FitModJ <- 'LJ =~ J1 + J2 + J3
            LJ ~ f*F + g*G + h*H
            F ~ b*D
            G ~ c*D + d*E
            H ~ e*E
            D ~~ a*E'
fitJ <- sem(FitModJ, sample.cov = covJ, sample.nobs = 100000)
summary(fitJ, std = TRUE)

FitModF <- 'LF =~ F1 + F2 + F3
            J ~ f*LF + g*G + h*H
            LF ~ b*D
            G ~ c*D + d*E
            H ~ e*E
            D ~~ a*E'
fitF <- sem(FitModF, sample.cov = covF, sample.nobs = 100000)
summary(fitF, std = TRUE)


FitModG <- 'LG =~ G1 + G2 + G3
            J ~ f*F + g*LG + h*H
            F ~ b*D
            LG ~ c*D + d*E
            H ~ e*E
            D ~~ a*E'
fitG <- sem(FitModG, sample.cov = covG, sample.nobs = 100000)
summary(fitG, std = TRUE)

FitModH <- 'LH =~ H1 + H2 + H3
            J ~ f*F + g*G + h*LH
            F ~ b*D
            G ~ c*D + d*E
            LH ~ e*E
            D ~~ a*E'
fitH <- sem(FitModH, sample.cov = covH, sample.nobs = 100000)
summary(fitH, std = TRUE)

FitModD <- 'LD =~ D1 + D2 + D3
            J ~ f*F + g*G + h*H
            F ~ b*LD
            G ~ c*LD + d*E
            H ~ e*E
            LD ~~ a*E'
fitD <- sem(FitModD, sample.cov = covD, sample.nobs = 100000)
summary(fitD, std = TRUE)

FitModE <- 'LE =~ E1 + E2 + E3
            J ~ f*F + g*G + h*H
            F ~ b*D
            G ~ c*D + d*LE
            H ~ e*LE
            D ~~ a*LE'
fitE <- sem(FitModE, sample.cov = covE, sample.nobs = 100000)
summary(fitE, std = TRUE)
################################################################################

# examine bias in explicit paths #################################################
summary(fitTrue, standardized = TRUE)
summary(fitJ, standardized = TRUE)
summary(fitF, standardized = TRUE)
summary(fitG, standardized = TRUE)
summary(fitH, standardized = TRUE)
summary(fitD, standardized = TRUE)
summary(fitE, standardized = TRUE)

#population RMSEA
sqrt((inspect(fitJ, "fit")["fmin"]*2)/inspect(fitJ,"fit")["df"])
sqrt((inspect(fitF, "fit")["fmin"]*2)/inspect(fitF,"fit")["df"])
sqrt((inspect(fitG, "fit")["fmin"]*2)/inspect(fitG,"fit")["df"])
sqrt((inspect(fitH, "fit")["fmin"]*2)/inspect(fitH,"fit")["df"])
sqrt((inspect(fitD, "fit")["fmin"]*2)/inspect(fitD,"fit")["df"])
sqrt((inspect(fitE, "fit")["fmin"]*2)/inspect(fitE,"fit")["df"])
################################################################################

# replicate C&P results ########################################################

# create measurement-error-laden versions of each variable: 
# to make this comparable to above: suppose 3 components per scale, and cov
# among components of .5. This implies reliability = .5*9 / 6 = .75

cov2 <- diag(7)
cov2[1:6,1:6] <- cov1
cov2[7,7] <- 1/3

sumMat <- matrix(0, 7, 13)
sumMat[1:7, 1:7] <- diag(7) #keep original 7 variables the same
sumMat[1:6, 8:13] <- diag(6) 
sumMat[7,8:13] <- 1 

cov3 <- t(sumMat) %*% cov2 %*% sumMat
colnames(cov3) <- c(colnames(cov1), "Err", paste(colnames(cov1), "E", sep = ""))

#run models with error-laden variables
ModErrJ <- 'JE ~ f*F + g*G + h*H
            F ~ b*D
            G ~ c*D + d*E
            H ~ e*E
            D ~~ a*E'

ModErrF <- 'J ~ f*FE + g*G + h*H
            FE ~ b*D
            G ~ c*D + d*E
            H ~ e*E
            D ~~ a*E'

ModErrG <- 'J ~ f*F + g*GE + h*H
            F ~ b*D
            GE ~ c*D + d*E
            H ~ e*E
            D ~~ a*E'

ModErrH <- 'J ~ f*F + g*G + h*HE
            F ~ b*D
            G ~ c*D + d*E
            HE ~ e*E
            D ~~ a*E'

ModErrD <- 'J ~ f*F + g*G + h*H
            F ~ b*DE
            G ~ c*DE + d*E
            H ~ e*E
            DE ~~ a*E'

ModErrE <- 'J ~ f*F + g*G + h*H
            F ~ b*D
            G ~ c*D + d*EE
            H ~ e*EE
            D ~~ a*EE'


fitErrJ <- sem(ModErrJ, sample.cov = cov3, sample.nobs = 100000)
fitErrF <- sem(ModErrF, sample.cov = cov3, sample.nobs = 100000)
fitErrG <- sem(ModErrG, sample.cov = cov3, sample.nobs = 100000)
fitErrH <- sem(ModErrH, sample.cov = cov3, sample.nobs = 100000)
fitErrD <- sem(ModErrD, sample.cov = cov3, sample.nobs = 100000)
fitErrE <- sem(ModErrE, sample.cov = cov3, sample.nobs = 100000)

#use summary function to get correct std.all estimates: 
summary(fitTrue, standardized = TRUE)
summary(fitErrJ, standardized = TRUE)
summary(fitErrF, standardized = TRUE)
summary(fitErrG, standardized = TRUE)
summary(fitErrH, standardized = TRUE)
summary(fitErrD, standardized = TRUE)
summary(fitErrE, standardized = TRUE)

#population RMSEA
sqrt((inspect(fitErrJ, "fit")["fmin"]*2)/inspect(fitErrJ,"fit")["df"])
sqrt((inspect(fitErrF, "fit")["fmin"]*2)/inspect(fitErrF,"fit")["df"])
sqrt((inspect(fitErrG, "fit")["fmin"]*2)/inspect(fitErrG,"fit")["df"])
sqrt((inspect(fitErrH, "fit")["fmin"]*2)/inspect(fitErrH,"fit")["df"])
sqrt((inspect(fitErrD, "fit")["fmin"]*2)/inspect(fitErrD,"fit")["df"])
sqrt((inspect(fitErrE, "fit")["fmin"]*2)/inspect(fitErrE,"fit")["df"])

#saturated path model: 
SatMod <- 'J ~ f*F + g*G + h*H + bJE*E + bJD*D
           F ~ b*D + bFE*E
           G ~ c*D + d*E
           H ~ e*E + bHD*D
           D ~~ a*E
           F ~~ sFG*G + sFH*H
           G ~~ sGH*H'

SatModErrJ <- 'JE ~ f*F + g*G + h*H + bJE*E + bJD*D
               F ~ b*D + bFE*E
               G ~ c*D + d*E
               H ~ e*E + bHD*D
               D ~~ a*E
               F ~~ sFG*G + sFH*H
               G ~~ sGH*H'

SatModErrF <- 'J ~ f*FE + g*G + h*H + bJE*E + bJD*D
               FE ~ b*D + bFE*E
               G ~ c*D + d*E
               H ~ e*E + bHD*D
               D ~~ a*E
               FE ~~ sFG*G + sFH*H
               G ~~ sGH*H'

SatModErrG <- 'J ~ f*F + g*GE + h*H + bJE*E + bJD*D
               F ~ b*D + bFE*E
               GE ~ c*D + d*E
               H ~ e*E + bHD*D
               D ~~ a*E
               F ~~ sFG*GE + sFH*H
               GE ~~ sGH*H'

SatModErrH <- 'J ~ f*F + g*G + h*HE + bJE*E + bJD*D
               F ~ b*D + bFE*E
               G ~ c*D + d*E
               HE ~ e*E + bHD*D
               D ~~ a*E
               F ~~ sFG*G + sFH*HE
               G ~~ sGH*HE'

SatModErrD <- 'J ~ f*F + g*G + h*H + bJE*E + bJD*DE
               F ~ b*DE + bFE*E
               G ~ c*DE + d*E
               H ~ e*E + bHD*DE
               DE ~~ a*E
               F ~~ sFG*G + sFH*H
               G ~~ sGH*H'

SatModErrE <- 'J ~ f*F + g*G + h*H + bJE*EE + bJD*D
               F ~ b*D + bFE*EE
               G ~ c*D + d*EE
               H ~ e*EE + bHD*D
               D ~~ a*EE
               F ~~ sFG*G + sFH*H
               G ~~ sGH*H'

fitSatErrJ <- sem(SatModErrJ, sample.cov = cov3, sample.nobs = 100000)
fitSatErrF <- sem(SatModErrF, sample.cov = cov3, sample.nobs = 100000)
fitSatErrG <- sem(SatModErrG, sample.cov = cov3, sample.nobs = 100000)
fitSatErrH <- sem(SatModErrH, sample.cov = cov3, sample.nobs = 100000)
fitSatErrD <- sem(SatModErrD, sample.cov = cov3, sample.nobs = 100000)
fitSatErrE <- sem(SatModErrE, sample.cov = cov3, sample.nobs = 100000)

#use summary function to get correct std.all estimates: 
summary(fitTrue, standardized = TRUE)
summary(fitSatErrJ, standardized = TRUE)
summary(fitSatErrF, standardized = TRUE)
summary(fitSatErrG, standardized = TRUE)
summary(fitSatErrH, standardized = TRUE)
summary(fitSatErrD, standardized = TRUE)
summary(fitSatErrE, standardized = TRUE)
 #results match those of C&P. 

################################################################################

# examine bias in implicit paths by adding implicit paths ###################### 
SatModLJ <- 'LJ =~ J1 + J2 + J3
            LJ ~ f*F + g*G + h*H + bJE*E + bJD*D
            F ~ b*D + bFE*E
            G ~ c*D + d*E
            H ~ e*E + bHD*D
            D ~~ a*E
            F ~~ sFG*G + sFH*H
            G ~~ sGH*H'

SatModLF <- 'LF =~ F1 + F2 + F3
            J ~ f*LF + g*G + h*H + bJE*E + bJD*D
            LF ~ b*D + bFE*E
            G ~ c*D + d*E
            H ~ e*E + bHD*D
            D ~~ a*E
            LF ~~ sFG*G + sFH*H
            G ~~ sGH*H'

SatModLG <- 'LG =~ G1 + G2 + G3
            J ~ f*F + g*LG + h*H + bJE*E + bJD*D
            F ~ b*D + bFE*E
            LG ~ c*D + d*E
            H ~ e*E + bHD*D
            D ~~ a*E
            F ~~ sFG*LG + sFH*H
            LG ~~ sGH*H'


SatModLH <- 'LH =~ H1 + H2 + H3
            J ~ f*F + g*G + h*LH + bJE*E + bJD*D
            F ~ b*D + bFE*E
            G ~ c*D + d*E
            LH ~ e*E + bHD*D
            D ~~ a*E
            F ~~ sFG*G + sFH*LH
            G ~~ sGH*LH'


SatModLD <- 'LD =~ D1 + D2 + D3
            J ~ f*F + g*G + h*H + bJE*E + bJD*LD
            F ~ b*LD + bFE*E
            G ~ c*LD + d*E
            H ~ e*E + bHD*LD
            LD ~~ a*E
            F ~~ sFG*G + sFH*H
            G ~~ sGH*H'


SatModLE <- 'LE =~ E1 + E2 + E3
            J ~ f*F + g*G + h*H + bJE*LE + bJD*D
            F ~ b*D + bFE*LE
            G ~ c*D + d*LE
            H ~ e*LE + bHD*D
            D ~~ a*LE
            F ~~ sFG*G + sFH*H
            G ~~ sGH*H'


fitSatLJ <- sem(SatModLJ, sample.cov = covJ, sample.nobs = 100000)
fitSatLF <- sem(SatModLF, sample.cov = covF, sample.nobs = 100000)
fitSatLG <- sem(SatModLG, sample.cov = covG, sample.nobs = 100000)
fitSatLH <- sem(SatModLH, sample.cov = covH, sample.nobs = 100000)
fitSatLD <- sem(SatModLD, sample.cov = covD, sample.nobs = 100000)
fitSatLE <- sem(SatModLE, sample.cov = covE, sample.nobs = 100000)

#use summary function to get correct std.all estimates: 
summary(fitSatLJ, standardized = TRUE)
summary(fitSatLF, standardized = TRUE)
summary(fitSatLG, standardized = TRUE)
summary(fitSatLH, standardized = TRUE)
summary(fitSatLD, standardized = TRUE)
summary(fitSatLE, standardized = TRUE)

################################################################################