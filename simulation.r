################################################################################
# Investigating Measurement Model Misspecification in Mediation Structures
# Author: Leonhard Volz
# Email: l.volz@uva.nl
# Date: 2021-09-03
################################################################################

# 0. libraries -----------------------------------------------------------------
library(lavaan)
library(R.utils)

# 1. Generating Functions ------------------------------------------------------
gen_adjacency_matrix <- function(size = 6,
                                 density = 1,
                                 entry = "single",
                                 shape = NA,
                                 layer_connected = FALSE) {
    if (!is.na(shape)) {
        if (sum(shape) != size) stop("Shape must compatible with square")
    }
    if (size && is.na(shape)) {
        if (entry == "single") {
            if (size == 6) {
                shape <- c(1, 2, 2, 1)
            } else if (size == 16) {
                shape <- c(1, 4, 6, 4, 1)
            } else if (size == 32) {
                shape <- c(1, 6, 6, 6, 6, 6, 1)
            }
        } else if (entry == "full") {
            if (size == 6) {
                shape <- c(3, 3)
            } else if (size == 16) {
                shape <- c(4, 4, 4, 4)
            } else if (size == 32) {
                shape <- c(6, 7, 6, 7, 6)
            }
        }
    }

    mat <- matrix(0, size + 2, size + 2)

    lower_layer_index <- 0
    lower_layer <- 1
    layer_index <- 1

    for (layer in shape) {
        for (j in (lower_layer_index + 1):(lower_layer_index + lower_layer)) {
            for (i in (layer_index + 1):(layer_index + layer)) {
                if (runif(1) < density) {
                    mat[i, j] <- 1
                }
            }
        }
        if (layer_connected) {
            for (j in (layer_index + 1):(layer_index + layer)) {
                for (i in (layer_index + 1):(layer_index + layer)) {
                    if (i != j && runif(1) < density) {
                        mat[i, j] <- 1
                    }
                }
            }
        }

        lower_layer_index <- layer_index
        lower_layer <- layer
        layer_index <- layer_index + layer
    }

    for (i in (lower_layer_index + 1):layer_index) {
        if (runif(1) < density) {
            mat[size + 2, i] <- 1
        }
    }

    rownames(mat) <- colnames(mat) <- c("X", paste0("m", 1:size), "Y")

    return(mat)
}


gen_loading_matrix <- function(mat,
                               x = .6,
                               m = .3,
                               y_prime = .6,
                               func = "constant") {
    if (length(x) == 1 && func == "constant") {
        mat[, 1] <- x * mat[, 1]
    }
    if (length(m) == 1 && func == "constant") {
        mat[, 2:(ncol(mat) - 1)] <- m * mat[2:(ncol(mat) - 1), ]
    }
    if (length(y_prime) == 1 && func == "constant") {
        mat[ncol(mat), ] <- y_prime * mat[ncol(mat), ]
    }

    return(mat)
}




generateCov.net.r <- function(a = .4, b = .4, beta = c(.3, .3, .3, .3), r = .8) {
    # total effect of X -> m4 = a:
    # r = reliability of observed indicators
    a. <- a / (beta[1] * beta[3] + beta[2] * beta[4])

    B <- matrix(0, 6, 6)
    B[2, 1] <- a.
    B[3, 2] <- beta[1]
    B[4, 2] <- beta[2]
    B[5, 3] <- beta[3]
    B[5, 4] <- beta[4]
    B[6, 5] <- b

    I <- diag(6)
    evec <- diag(B %*% t(B))
    evec[5] <- evec[5] + 2 * beta[1] * beta[2] * beta[3] * beta[4]

    Psi <- diag(6) - diag(evec)
    Sigma <- solve(I - B) %*% Psi %*% t(solve(I - B))
    # order = x, m1, m2, m3, m4, y

    # add a column for M (sum of m1-m4)
    sumMat <- matrix(0, 6, 7)
    sumMat[1:6, 1:6] <- diag(6)
    sumMat[2:5, 7] <- 1

    Sigma <- t(sumMat) %*% Sigma %*% sumMat
    # order = x, m1, m2, m3, m4, y, M

    # create unreliable indicators of m1-m4, call them y1-y4:
    Lambda <- matrix(0, 11, 7)
    diag(Lambda[1:7, 1:7]) <- 1
    diag(Lambda[8:11, 1:4]) <- sqrt(r)
    Theta <- matrix(0, 11, 11)
    diag(Theta[8:11, 8:11]) <- 1 - r

    Sigma2 <- Lambda %*% Sigma %*% t(Lambda) + Theta
    # order = x, m1, m2, m3, m4, y, M, y1, y2, y3, y4

    # add a column for M (sum of y1-y4)
    sumMat <- matrix(0, 11, 12)
    sumMat[1:11, 1:11] <- diag(11)
    sumMat[8:11, 12] <- 1

    Sigma3 <- t(sumMat) %*% Sigma2 %*% sumMat
    # order = x, m1, m2, m3, m4, y, M, y1, y2, y3, y4, M2

    colnames(Sigma3) <- c("X", "m1", "m2", "m3", "m4", "Y", "M", "y1", "y2", "y3", "y4", "M2")

    return(Sigma3)
}
