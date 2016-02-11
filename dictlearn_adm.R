# ADM FOR SPARSE DICTIONARY LEARNING ---------------------------
# dictlearn_adm.R
#
# Tue Feb 9 2016
# From Sun, Qu, and Wright, "Complete Dictionary Recovery over the Sphere", ICML 2015
#
# Code based off: https://github.com/sunju/dl_focm/tree/master/image_experiment_focm
#

## Read in image data (probably better to download images locally for testing, but URL should also work)
library(pixmap)
im <- read.pnm("https://github.com/sunju/dl_focm/blob/master/image_experiment_focm/test_slate/25.pgm") # can be any number between 25-36

image_to_patches <- function(dat, patch_size){
    # patch_size <- 8
    dim.out <- dim(dat) / patch_size
    outY <- matrix(0, nrow = patch_size*patch_size, ncol = prod(dim.out))

    k <- 1
    for(i in 1:dim.out[1]){
        for(j in 1:dim.out[2]){
            outY[, k] <- as.vector(dat[((i-1)*patch_size + 1):(i*patch_size), ((j-1)*patch_size + 1):(j*patch_size)])
            k <- k+1
        }
    }

    outY
}

proj_ortho_group <- function(m){
    m.svd <- svd(m)
    m.svd$u %*% t(m.svd$v)
}

prox_l1 <- function(m, lambda){
    sign(m) * pmax(m - lambda, 0)
}

run_adm_orthobasis <- function(Y, initA,
                               lambda, tol, iters,
                               verbose = TRUE, iter.out = 100){
    oldA <- initA

    frob_diffs <- rep(NA, iters)
    obj_vals <- rep(NA, iters)
    for(iter in 1:iters){

        ## Show and plot progress
        if(verbose && iter%%iter.out == 0){
            cat(sprintf("Working on %d/%d...\n", iter, iters))
            par(mfrow=c(2,1), mai=c(1.1,1.1,0.1,0.1))
            plot(frob_diffs, pch=16, cex=0.5, xlab="")
            plot(obj_vals, pch=16, cex=0.5, col="blue")
        }

        ## Update (A, X)
        X <- prox_l1(t(oldA) %*% Y, lambda)
        A <- proj_ortho_group(Y %*% t(X))

        ## Compute difference between estimates of A
        frob_diffs[iter] <- sum((A - oldA)^2)
        if(frob_diffs[iter] < tol) break

        ## Compute value of objective function
        reconstruction <- A %*% X
        obj_vals[iter] <- 0.5 * sum((reconstruction - Y)^2) + lambda * sum(abs(X))

        oldA <- A
    }

    list(A = A, X = X, diffA = frob_diffs, obj = obj_vals, lambda = lambda)
}

#----------------------------------------------------------------------------------------------------
# RUN THE MAIN ALGORITHM
#----------------------------------------------------------------------------------------------------

## Data munging
dat <- getChannels(im)          # convert pixmap to matrix
Y <- image_to_patches(dat, 8)   # transform matrix to "image patches" (Section 1.2)
nn <- nrow(Y)                   # dimension of patches
pp <- ncol(Y)                   # number of samples
mm <- nn                        # size of dictionary (complete dictionary => m=n)

## Set up initial parameters
iters <- 1000   # max number of iterations
lambda <- .2    # tau in original code
tol <- 1e-5     # tolerance parameter for convergence
initA <- proj_ortho_group(matrix(rnorm(nn*mm), nrow = nn))  # generate a random orthobasis

## Run the ADM algorithm
adm.out <- run_adm_orthobasis(Y, initA, lambda, tol, iters)

## Plot the errors and the objective
plot(adm.out$diffA, pch=16, cex=0.3)
plot(adm.out$obj, pch=16, cex=0.5, col="blue")

## How well does the reconstructed image AX match the input Y?
reconstruction <- adm.out$A %*% adm.out$X
summary(as.vector(abs(reconstruction - Y)))
