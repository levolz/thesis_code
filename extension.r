################################################################################
# R code to extend Rhemtulla et al. (2021)
# Author: Leonhard Volz
# Email: l.volz@uva.nl
# Date: 2021-09-03
################################################################################

# 0. libraries -----------------------------------------------------------------
library(lavaan)
library(R.utils)

# 1. Models --------------------------------------------------------------------
## LV model
model_base <- "
    Mb =~ a * X
    Y ~ b * Mb
    ab := a*b
"

model_lv <- "
    M =~ m1 + m2 + m3 + m4 + m5 + m6
    M ~ a * X
    Y ~ b * M
    ab := a*b
"

model_full <- "
    M =~ m1 + m2 + m3 + m4 + m5 + m6
    m1 ~ X
    m2 ~ X
    m3 ~ X
    m4 ~ X
    m5 ~ X
    m6 ~ X
"

# 2. Generate Functions ---------------------------------------------------------
## Generate Covariance Matrix for Composite/Causal Indicators
generatecov_ci <- function(
    a = .6,
    b = .4,
    avec = c(.35, .4, .45, .55, .6, .65),
    # covm = NA,
    res_m = 0,
    r = NA) {
    amat <- avec %*% t(avec)

    gamma <- a / sum(avec)
    covm <- (1 - 6 * gamma^2 - res_m) / (30 * gamma^2)

    evec <- 1 - avec^2
    emat <- covm - amat

    beta <- matrix(0, 9, 9)
    beta[1:6, 9] <- avec
    beta[7, 1:6] <- gamma
    beta[8, 7] <- b
    I <- diag(9)
    psi <- matrix(0, 9, 9)
    psi[1:6, 1:6] <- emat
    diag(psi) <- c(evec, res_m, (1 - b^2), 1)

    sigma <- solve(I - beta) %*% psi %*% t(solve(I - beta))

    if (is.na(r)) {
        return(sigma)
    }

    lambda <- matrix(0, 15, 9)
    diag(lambda[1:9, 1:9]) <- 1
    diag(lambda[10:15, 1:6]) <- sqrt(r)
    theta <- matrix(0, 15, 15)
    diag(theta[10:15, 10:15]) <- 1 - r

    sigma_i <- lambda %*% sigma %*% t(lambda) + theta

    sumMat <- matrix(0, 15, 16)
    sumMat[1:15, 1:15] <- diag(15)
    sumMat[10:15, 16] <- 1

    sigma_ii <- t(sumMat) %*% sigma_i %*% sumMat

    colnames(sigma_ii) <- c(paste0("m", 1:6), "M", "Y", "X", paste0("y", 1:6), "Mb")
    rownames(sigma_ii) <- c(paste0("m", 1:6), "M", "Y", "X", paste0("y", 1:6), "Mb")

    return(sigma_ii)
}

### Example
Sigma <- generatecov_ci(r = .8)
fit_lv <- sem(model = model_lv, sample.cov = Sigma, sample.nobs = 10000, std.lv = TRUE)
fit_full <- sem(model = model_full, sample.cov = Sigma, sample.nobs = 10000, std.lv = TRUE)

summary(fit_lv)
summary(fit_full)

## Generate Covariance Matrix for Latent Variable Model
# generate_cov_lv <- function(
#     a = .6,
#     b = .4,
#     avec = c(.35, .4, .45, .55, .6, .65),
#     res_m = 0,
#     r = NA) {
#     gamma <- a / sum(avec)
#     covm <- (1 - 6 * gamma^2 - res_m) / (30 * gamma^2)
#
#     evec <- 1 - avec^2
#     emat <- covm - avec %*% t(avec)
#
#     beta <- matrix(0, 7, 7)
#     beta[1:6, 7] <- avec
#     beta[7, 1:6] <- gamma
#     beta[7, 7] <- b
#     I <- diag(7)
#     psi <- matrix(0, 7, 7)
#     psi[1:6, 1:6] <- emat
#     diag(psi) <- c(evec, res_m, (1 - b^2))
#
#     sigma <- solve(I - beta) %*% psi %*% t(solve(I - beta))
#
#     if (!r) {
#         return(sigma)
#     }
#
#     lambda <- matrix(0, 13, 7)
#     diag(lambda[1:7, 1:7]) <- 1
#     diag(lambda[8:13, 1:6]) <- sqrt(r)
#     theta <- matrix(0, 13, 13)
#     diag(theta[8:13, 8:13]) <- 1 - r
#
#     sigma_i <- lambda %*% sigma %*% t(lambda) + theta
#
#     sumMat <- matrix(0, 13, 14)
#     sumMat[1:13, 1:13] <- diag(13)
#     sumMat[8:13, 14] <- 1
#
#     sigma_ii <- t(sumMat) %*% sigma_i %*% sumMat
#
#     colnames(sigma_ii) <- c(paste0("m", 1:6), "M", "Y", "X", paste0("y", 1:6))
#     rownames(sigma_ii) <- c(paste0("m", 1:6), "M", "Y", "X", paste0("y", 1:6))
#
#     return(sigma_ii)
# }

## Generate MIMIC Model
generate_cov_mimic <- function(
    a = .6,
    b = .4,
    avec = c(.5, .4, .3),
    lamvec = c(.5, .6, .7),
    covm = NA,
    res.M = 0,
    r = .8) {
    amat <- avec %*% t(avec)

    if (is.na(covm)) {
        gamma <- a / sum(avec)
        covm <- (1 - 3 * gamma^2 - res.M) / (6 * gamma^2)
    } else {
        gamma <- 1 / (3 * covm + 3)
        a <- gamma * sum(avec)
    }

    evec <- 1 - avec^2
    emat <- covm - amat

    beta <- matrix(0, 9, 9)
    beta[1:3, 9] <- avec
    beta[7, 1:3] <- gamma
    beta[4:6, 7] <- lamvec
    beta[8, 7] <- b

    I <- diag(9)
    psi <- matrix(0, 9, 9)
    psi[1:3, 1:3] <- emat
    diag(psi) <- c(evec, 1 - lamvec^2, res.M, (1 - b^2), 1)

    sigma <- solve(I - beta) %*% psi %*% t(solve(I - beta))

    if (is.na(r)) {
        return(sigma)
    }

    lambda <- matrix(0, 12, 9)
    diag(lambda[1:9, 1:9]) <- 1
    diag(lambda[10:12, 1:3]) <- sqrt(r)
    theta <- matrix(0, 12, 12)
    diag(theta[10:12, 10:12]) <- 1 - r

    sigma_i <- lambda %*% sigma %*% t(lambda) + theta

    sumMat <- matrix(0, 12, 13)
    sumMat[1:12, 1:12] <- diag(12)
    sumMat[10:12, 13] <- 1

    sigma_ii <- t(sumMat) %*% sigma_i %*% sumMat

    colnames(sigma_ii) <- c(paste0("m", 1:6), "M", "Y", "X", paste0("y", 1:3), "Mb")
    rownames(sigma_ii) <- c(paste0("m", 1:6), "M", "Y", "X", paste0("y", 1:3), "Mb")

    return(sigma_ii)
}


## Generate Covariance Matrix for Network Model
# generate_cov_net <- function(
#     a = .6,
#     b = .4,
#     beta = matrix(c())
